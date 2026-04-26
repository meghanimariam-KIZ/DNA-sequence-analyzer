from __future__ import annotations

from pathlib import Path
from typing import List

from Bio import SeqIO
from Bio.Align import PairwiseAligner


SAMPLE_DB_PATH = Path(__file__).resolve().parents[2] / "sample_data" / "sample_database.fasta"


def _parse_sample_description(description: str) -> tuple[str, str]:
    parts = description.split(maxsplit=1)
    accession = parts[0]
    organism = parts[1] if len(parts) > 1 else accession
    return accession, organism


def _load_demo_database() -> list[dict]:
    records = []
    for record in SeqIO.parse(SAMPLE_DB_PATH, "fasta"):
        accession, organism = _parse_sample_description(record.description)
        records.append(
            {
                "accession_id": accession,
                "organism_name": organism,
                "sequence": str(record.seq).upper(),
            }
        )
    return records


def _alignment_metrics(query: str, subject: str) -> dict:
    """Compute demo BLAST-like metrics with pairwise global alignment.

    This is intentionally simple and local. It is not a replacement for NCBI
    BLAST, but it gives students realistic fields to learn from.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(query, subject)[0]
    aligned_query = alignment.aligned[0]
    aligned_subject = alignment.aligned[1]

    matches = 0
    aligned_bases = 0
    query_covered = set()

    for q_block, s_block in zip(aligned_query, aligned_subject):
        q_start, q_end = q_block
        s_start, s_end = s_block
        block_length = min(q_end - q_start, s_end - s_start)
        for offset in range(block_length):
            q_index = q_start + offset
            s_index = s_start + offset
            query_covered.add(q_index)
            if query[q_index] == subject[s_index]:
                matches += 1
            aligned_bases += 1

    percent_identity = round((matches / aligned_bases * 100) if aligned_bases else 0, 2)
    query_coverage = round((len(query_covered) / len(query) * 100) if query else 0, 2)
    score = round(float(alignment.score), 2)
    e_value = "{:.2e}".format(max(1e-180, 10 ** (-(score / 25))))

    return {
        "percent_identity": percent_identity,
        "query_coverage": query_coverage,
        "alignment_score": score,
        "e_value": e_value,
    }


def runBlastSearch(sequence: str) -> List[dict]:
    """Run similarity search.

    Current mode: local demo database.

    Future extension:
    Replace the body of this function with a call to NCBI BLAST+ or NCBI's
    remote BLAST service. The API response can be mapped to the same output
    keys so the frontend does not need to change.
    """
    results = []
    for record in _load_demo_database():
        metrics = _alignment_metrics(sequence, record["sequence"])
        results.append(
            {
                "organism_name": record["organism_name"],
                "accession_id": record["accession_id"],
                "sequence": record["sequence"],
                **metrics,
            }
        )

    return sorted(
        results,
        key=lambda item: (item["percent_identity"], item["query_coverage"], item["alignment_score"]),
        reverse=True,
    )[:6]

