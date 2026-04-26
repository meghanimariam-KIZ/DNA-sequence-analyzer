from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from io import StringIO
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


DNA_ALPHABET = set("ACGTN")
CODON_TABLE = 1
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class CleanSequence:
    sequence: str
    header: str | None
    removed_characters: int


def parse_and_clean_sequence(raw_text: str) -> CleanSequence:
    """Parse FASTA/plain text and keep only DNA characters.

    Beginner note:
    DNA sequence files often contain headers, spaces, line numbers, or accidental
    punctuation. This function extracts sequence letters and removes anything
    that is not A, C, G, T, or N.
    """
    if not raw_text or not raw_text.strip():
        raise ValueError("Please paste or upload a DNA sequence.")

    text = raw_text.strip()
    header = None
    sequence_parts: List[str] = []

    if text.startswith(">"):
        records = list(SeqIO.parse(StringIO(text), "fasta"))
        if not records:
            raise ValueError("The FASTA input could not be parsed.")
        record = records[0]
        header = record.description
        raw_sequence = str(record.seq)
    else:
        raw_sequence = text

    raw_sequence = raw_sequence.upper()
    for character in raw_sequence:
        if character in DNA_ALPHABET:
            sequence_parts.append(character)

    sequence = "".join(sequence_parts)
    removed = len(raw_sequence) - len(sequence)

    if not sequence:
        raise ValueError("No valid DNA bases were found. Use A, C, G, T, or N.")
    if len(sequence) < 9:
        raise ValueError("Sequence is too short for useful ORF and similarity analysis.")

    return CleanSequence(sequence=sequence, header=header, removed_characters=max(removed, 0))


def gc_at_content(sequence: str) -> Tuple[float, float]:
    counted = [base for base in sequence if base in {"A", "C", "G", "T"}]
    if not counted:
        return 0.0, 0.0
    gc = counted.count("G") + counted.count("C")
    at = counted.count("A") + counted.count("T")
    return round(gc / len(counted) * 100, 2), round(at / len(counted) * 100, 2)


def estimate_molecular_weight(sequence: str) -> float:
    """Estimate single-stranded DNA molecular weight in Daltons."""
    weights = {"A": 313.21, "T": 304.2, "C": 289.18, "G": 329.21, "N": 308.95}
    return round(sum(weights.get(base, 0.0) for base in sequence), 2)


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def transcribe(sequence: str) -> str:
    return str(Seq(sequence).transcribe())


def translate_sequence(sequence: str, frame: int = 0) -> str:
    trimmed = sequence[frame:]
    usable_length = len(trimmed) - (len(trimmed) % 3)
    if usable_length <= 0:
        return ""
    return str(Seq(trimmed[:usable_length]).translate(table=CODON_TABLE, to_stop=False))


def codon_usage(sequence: str) -> Dict[str, int]:
    counts: Counter[str] = Counter()
    usable_length = len(sequence) - (len(sequence) % 3)
    for index in range(0, usable_length, 3):
        codon = sequence[index:index + 3]
        if len(codon) == 3:
            counts[codon] += 1
    return dict(sorted(counts.items()))


def find_orfs_in_frame(sequence: str, frame_number: int, strand: str) -> List[dict]:
    """Find ORFs that start with ATG and end at the first in-frame stop codon."""
    orfs: List[dict] = []
    frame_offset = abs(frame_number) - 1
    sequence_length = len(sequence)

    for start in range(frame_offset, sequence_length - 2, 3):
        codon = sequence[start:start + 3]
        if codon != "ATG":
            continue

        for stop in range(start + 3, sequence_length - 2, 3):
            stop_codon = sequence[stop:stop + 3]
            if stop_codon in STOP_CODONS:
                dna = sequence[start:stop + 3]
                protein = str(Seq(dna).translate(table=CODON_TABLE, to_stop=True))
                orfs.append(
                    {
                        "frame": frame_number,
                        "strand": strand,
                        "start": start + 1,
                        "end": stop + 3,
                        "length": len(dna),
                        "start_codon": "ATG",
                        "stop_codon": stop_codon,
                        "dna_sequence": dna,
                        "protein_sequence": protein,
                    }
                )
                break

    return orfs


def find_orfs_all_frames(sequence: str) -> Tuple[List[dict], List[dict], dict | None]:
    forward_orfs: List[dict] = []
    reverse_orfs: List[dict] = []
    rc_sequence = reverse_complement(sequence)

    for frame in (1, 2, 3):
        forward_orfs.extend(find_orfs_in_frame(sequence, frame, "+"))
        reverse_orfs.extend(find_orfs_in_frame(rc_sequence, -frame, "-"))

    all_orfs = sorted(forward_orfs + reverse_orfs, key=lambda item: item["length"], reverse=True)
    longest = all_orfs[0] if all_orfs else None

    frame_info = []
    for frame in (1, 2, 3):
        translated = translate_sequence(sequence, frame - 1)
        frame_info.append(
            {
                "frame": frame,
                "strand": "+",
                "protein_preview": translated[:120],
                "protein_length": len(translated),
                "orf_count": len([orf for orf in forward_orfs if orf["frame"] == frame]),
            }
        )
    for frame in (1, 2, 3):
        translated = translate_sequence(rc_sequence, frame - 1)
        frame_info.append(
            {
                "frame": -frame,
                "strand": "-",
                "protein_preview": translated[:120],
                "protein_length": len(translated),
                "orf_count": len([orf for orf in reverse_orfs if orf["frame"] == -frame]),
            }
        )

    return all_orfs, frame_info, longest


def analyze_sequence(raw_text: str) -> dict:
    cleaned = parse_and_clean_sequence(raw_text)
    sequence = cleaned.sequence
    gc_percent, at_percent = gc_at_content(sequence)
    all_orfs, frame_info, longest_orf = find_orfs_all_frames(sequence)

    return {
        "input_header": cleaned.header,
        "cleaned_sequence": sequence,
        "removed_characters": cleaned.removed_characters,
        "length": len(sequence),
        "gc_percent": gc_percent,
        "at_percent": at_percent,
        "molecular_weight_da": estimate_molecular_weight(sequence),
        "reverse_complement": reverse_complement(sequence),
        "rna_sequence": transcribe(sequence),
        "protein_translation": translate_sequence(sequence),
        "codon_usage": codon_usage(sequence),
        "orfs": all_orfs[:100],
        "longest_orf": longest_orf,
        "reading_frames": frame_info,
    }

