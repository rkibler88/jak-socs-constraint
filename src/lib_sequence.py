"""Shared utilities for the JAK/SOCS inter-protein coevolution pipeline.

Covers sequence retrieval, alignment, gap-aware coordinate mapping between
canonical UniProt numbering and alignment columns, and effective-depth
statistics. Kept in one module so every pipeline step uses the same
conventions - in particular the same N_eff definition, since that number gates
the whole analysis.
"""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

LOG = logging.getLogger("dcalib")

#: Canonical amino-acid alphabet plus gap. Index 0 is the gap symbol, matching
#: the convention used by plmc and CCMpred.
ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {c: i for i, c in enumerate(ALPHABET)}
Q = len(ALPHABET)

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_ACCESSIONS = "https://rest.uniprot.org/uniprotkb/accessions"


# ---------------------------------------------------------------- FASTA I/O

def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a FASTA file into an ordered {header: sequence} mapping."""
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def write_fasta(records: dict[str, str], path: str | Path, width: int = 60) -> None:
    """Write {header: sequence} to FASTA with fixed line width."""
    with open(path, "w") as fh:
        for header, seq in records.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# ------------------------------------------------------- UniProt retrieval

def fetch_sequences(accessions: list[str], batch: int = 150,
                    contact_email: str | None = None) -> dict[str, str]:
    """Fetch canonical sequences for ``accessions`` from UniProt.

    Queries are batched as accession disjunctions against the stream endpoint,
    which is markedly faster than one request per accession and stays within
    UniProt's fair-use expectations.
    """
    headers = {"Accept": "text/plain"}
    if contact_email:
        headers["User-Agent"] = f"jak-socs-coevolution ({contact_email})"

    out: dict[str, str] = {}
    unique = list(dict.fromkeys(accessions))
    for start in range(0, len(unique), batch):
        chunk = unique[start:start + batch]
        # The dedicated bulk-accessions endpoint avoids the URL-length limit
        # that an accession disjunction hits at this batch size.
        params = {"accessions": ",".join(chunk), "format": "fasta"}
        url = f"{UNIPROT_ACCESSIONS}?{urllib.parse.urlencode(params)}"
        body = _get(url, headers)
        for header, seq in _parse_fasta_string(body).items():
            # UniProt FASTA headers look like sp|ACC|NAME description
            parts = header.split("|")
            acc = parts[1] if len(parts) > 2 else header.split()[0]
            out[acc] = seq.upper()
        LOG.info("fetched %d/%d sequences", len(out), len(unique))
    missing = set(unique) - set(out)
    if missing:
        LOG.warning("%d accessions returned no sequence (e.g. %s)",
                    len(missing), sorted(missing)[:5])
    return out


def _get(url: str, headers: dict[str, str], retries: int = 4) -> str:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read().decode("utf-8")
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as err:
            if isinstance(err, urllib.error.HTTPError) and 400 <= err.code < 500:
                raise
            last = err
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"request failed after {retries} attempts: {last}")


def _parse_fasta_string(text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header, chunks = line[1:], []
        elif line:
            chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks)
    return records


# ------------------------------------------------------------- Alignment

def run_mafft(records: dict[str, str], out_path: str | Path,
              mode: str = "auto", threads: int = 8) -> dict[str, str]:
    """Align ``records`` with MAFFT and return the aligned records.

    mode ``auto`` uses ``--auto``; ``einsi`` uses ``--genafpair --maxiterate
    1000``, appropriate for sequences with large gaps such as the SOCS set.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tmp:
        write_fasta(records, tmp.name)
        tmp_in = tmp.name

    flags = {"auto": ["--auto"],
             "einsi": ["--genafpair", "--maxiterate", "1000"],
             "nsi": ["--retree", "2", "--maxiterate", "1000"],
             "fast": ["--retree", "1", "--maxiterate", "0"]}[mode]
    cmd = ["mafft", *flags, "--anysymbol", "--thread", str(threads), tmp_in]
    LOG.info("running: %s", " ".join(cmd[:-1]) + f" <{len(records)} seqs>")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed ({proc.returncode}): {proc.stderr[-2000:]}")
    aligned = _parse_fasta_string(proc.stdout)
    aligned = {k: v.upper() for k, v in aligned.items()}
    write_fasta(aligned, out_path)
    lengths = {len(v) for v in aligned.values()}
    assert len(lengths) == 1, f"MAFFT returned ragged alignment: {sorted(lengths)[:5]}"
    return aligned


# -------------------------------------------------- Coordinate mapping

@dataclass(frozen=True)
class CoordinateMap:
    """Bijection between canonical residue numbers and alignment columns.

    ``human_to_msa`` maps 1-based canonical residue number to 0-based alignment
    column. ``msa_to_human`` is its inverse, defined only on columns where the
    reference is not gapped. Built from the aligned reference row, so
    insertions, deletions and gaps are handled by construction rather than by
    an offset.
    """

    reference: str
    human_to_msa: dict[int, int]
    msa_to_human: dict[int, int]
    n_columns: int

    def validate(self, canonical_seq: str, aligned_ref: str) -> None:
        """Round-trip every canonical residue and assert identity is preserved.

        Raises AssertionError on any mismatch. This is the check whose absence
        let a previous pipeline slice a 'SOCS' block that contained no SOCS
        residues.
        """
        assert len(self.human_to_msa) == len(canonical_seq), (
            f"{self.reference}: mapped {len(self.human_to_msa)} of "
            f"{len(canonical_seq)} canonical residues"
        )
        for pos in range(1, len(canonical_seq) + 1):
            col = self.human_to_msa[pos]
            assert aligned_ref[col] == canonical_seq[pos - 1], (
                f"{self.reference}: residue {pos} "
                f"({canonical_seq[pos-1]}) maps to column {col} "
                f"which holds {aligned_ref[col]!r}"
            )
            assert self.msa_to_human[col] == pos, (
                f"{self.reference}: round trip failed at residue {pos}"
            )


def build_coordinate_map(reference: str, aligned_ref: str,
                         canonical_seq: str | None = None) -> CoordinateMap:
    """Build and validate a CoordinateMap from an aligned reference row."""
    human_to_msa: dict[int, int] = {}
    msa_to_human: dict[int, int] = {}
    pos = 0
    for col, ch in enumerate(aligned_ref):
        if ch not in ("-", ".", " "):
            pos += 1
            human_to_msa[pos] = col
            msa_to_human[col] = pos
    cmap = CoordinateMap(reference, human_to_msa, msa_to_human, len(aligned_ref))
    if canonical_seq is not None:
        cmap.validate(canonical_seq, aligned_ref)
    return cmap


def project_range(cmap: CoordinateMap, start: int, end: int) -> list[int]:
    """Project an inclusive canonical residue range to alignment columns."""
    missing = [p for p in range(start, end + 1) if p not in cmap.human_to_msa]
    assert not missing, (
        f"{cmap.reference}: residues {missing[:5]} of range {start}-{end} "
        "are absent from the coordinate map"
    )
    return [cmap.human_to_msa[p] for p in range(start, end + 1)]


# ------------------------------------------------------- Effective depth

def encode_msa(records: dict[str, str]) -> np.ndarray:
    """Encode aligned records as an (N, L) int8 array over ALPHABET.

    Any symbol outside ALPHABET (X, B, Z, U, O, ambiguity codes) is mapped to
    the gap index, which is how plmc and CCMpred treat unknown residues.
    """
    seqs = list(records.values())
    lengths = {len(s) for s in seqs}
    assert len(lengths) == 1, f"ragged alignment: {sorted(lengths)[:5]}"
    n, length = len(seqs), lengths.pop()
    arr = np.zeros((n, length), dtype=np.int8)
    for i, seq in enumerate(seqs):
        arr[i] = [AA_INDEX.get(c, 0) for c in seq]
    return arr


def sequence_weights(msa: np.ndarray, theta: float = 0.8,
                     block: int = 512) -> np.ndarray:
    """Redundancy-downweighting weights at identity threshold ``theta``.

    Two sequences are neighbours when the fraction of alignment columns at
    which they carry the identical symbol is at least ``theta``, with gaps
    counted as a symbol. Each sequence's weight is the reciprocal of its
    neighbour count including itself. This is the standard DCA convention
    (Morcos et al. 2011) and the one plmc implements, so N_eff here is
    comparable to the value plmc reported for the previous run.
    """
    n, length = msa.shape
    counts = np.zeros(n, dtype=np.float64)
    for start in range(0, n, block):
        stop = min(start + block, n)
        # (stop-start, n) identity fractions
        same = (msa[start:stop, None, :] == msa[None, :, :]).sum(axis=2) / length
        counts[start:stop] = (same >= theta).sum(axis=1)
    return 1.0 / counts


def effective_depth(msa: np.ndarray, theta: float = 0.8) -> float:
    """N_eff: the sum of redundancy-downweighted sequence weights."""
    return float(sequence_weights(msa, theta).sum())


def gap_fraction(msa: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """Per-column gap fraction, optionally weighted by sequence weights."""
    is_gap = (msa == 0)
    if weights is None:
        return is_gap.mean(axis=0)
    w = weights / weights.sum()
    return (is_gap * w[:, None]).sum(axis=0)


@dataclass
class DepthReport:
    """Depth diagnostics for one paired alignment."""

    dataset: str
    n_sequences: int
    l_jak: int
    l_socs: int
    neff_80: float

    @property
    def total_length(self) -> int:
        return self.l_jak + self.l_socs

    @property
    def neff_over_l(self) -> float:
        return self.neff_80 / self.total_length

    @property
    def ovchinnikov_threshold(self) -> float:
        """N_eff needed for reliable inter-protein couplings: mean chain length.

        Ovchinnikov, Kamisetty & Baker, eLife 2014 (10.7554/eLife.02030).
        """
        return self.total_length / 2.0

    @property
    def passes_gate(self) -> bool:
        return self.neff_80 > self.ovchinnikov_threshold

    @property
    def gate_ratio(self) -> float:
        """N_eff as a multiple of the required threshold; >1 passes."""
        return self.neff_80 / self.ovchinnikov_threshold
