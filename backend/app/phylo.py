from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Tuple

from Bio import pairwise2
from Bio import Phylo
from Bio.Phylo.Newick import Clade, Tree


def simple_star_msa(sequences: Dict[str, str]) -> Dict[str, str]:
    """Create a beginner-friendly star multiple sequence alignment.

    The query sequence is used as the center. Each hit is globally aligned to
    the query, then gap positions are merged into one shared alignment. This is
    useful for teaching, but production phylogenetics should use tools such as
    MAFFT, MUSCLE, Clustal Omega, or validated BioPython/scikit-bio workflows.
    """
    names = list(sequences.keys())
    reference_name = names[0]
    reference = sequences[reference_name]

    parsed_alignments = {
        reference_name: {
            "bases": list(reference),
            "insertions": {index: [] for index in range(len(reference) + 1)},
        }
    }

    for name in names[1:]:
        alignment = pairwise2.align.globalms(reference, sequences[name], 2, -1, -4, -0.5, one_alignment_only=True)[0]
        aligned_reference, aligned_sequence = alignment.seqA, alignment.seqB
        bases = ["-"] * len(reference)
        insertions = {index: [] for index in range(len(reference) + 1)}
        reference_position = 0

        for ref_base, seq_base in zip(aligned_reference, aligned_sequence):
            if ref_base == "-":
                insertions[reference_position].append(seq_base)
            else:
                bases[reference_position] = seq_base
                reference_position += 1

        parsed_alignments[name] = {"bases": bases, "insertions": insertions}

    insertion_widths = {
        index: max(len(parsed_alignments[name]["insertions"][index]) for name in names)
        for index in range(len(reference) + 1)
    }

    msa: Dict[str, str] = {}
    for name in names:
        pieces: List[str] = []
        parsed = parsed_alignments[name]
        for index in range(len(reference)):
            inserted = parsed["insertions"][index]
            pieces.extend(inserted)
            pieces.extend("-" for _ in range(insertion_widths[index] - len(inserted)))
            pieces.append(parsed["bases"][index])

        inserted = parsed["insertions"][len(reference)]
        pieces.extend(inserted)
        pieces.extend("-" for _ in range(insertion_widths[len(reference)] - len(inserted)))
        msa[name] = "".join(pieces)

    return msa


def _p_distance(seq_a: str, seq_b: str) -> float:
    comparable = [(a, b) for a, b in zip(seq_a, seq_b) if a != "-" and b != "-"]
    if not comparable:
        return 1.0
    mismatches = sum(1 for a, b in comparable if a != b)
    return round(mismatches / len(comparable), 5)


def build_distance_matrix(aligned: Dict[str, str]) -> Dict[str, Dict[str, float]]:
    names = list(aligned.keys())
    matrix = {name: {other: 0.0 for other in names} for name in names}

    for name_a, name_b in combinations(names, 2):
        distance = _p_distance(aligned[name_a], aligned[name_b])
        matrix[name_a][name_b] = distance
        matrix[name_b][name_a] = distance

    return matrix


def _average_distance(cluster_a: List[str], cluster_b: List[str], original: Dict[str, Dict[str, float]]) -> float:
    distances = [original[a][b] for a in cluster_a for b in cluster_b]
    return sum(distances) / len(distances)


def build_upgma_tree(sequences: Dict[str, str]) -> Tuple[str, dict, Dict[str, Dict[str, float]], Dict[str, str]]:
    """Build a simple UPGMA tree from pairwise p-distances."""
    alignment = simple_star_msa(sequences)
    distance_matrix = build_distance_matrix(alignment)
    clusters: Dict[str, List[str]] = {name: [name] for name in sequences}
    clades: Dict[str, Clade] = {name: Clade(name=name) for name in sequences}
    heights: Dict[str, float] = {name: 0.0 for name in sequences}

    while len(clusters) > 1:
        names = list(clusters.keys())
        best_pair = None
        best_distance = float("inf")

        for name_a, name_b in combinations(names, 2):
            distance = _average_distance(clusters[name_a], clusters[name_b], distance_matrix)
            if distance < best_distance:
                best_pair = (name_a, name_b)
                best_distance = distance

        assert best_pair is not None
        name_a, name_b = best_pair
        new_name = f"({name_a},{name_b})"
        new_height = best_distance / 2

        clade_a = clades[name_a]
        clade_b = clades[name_b]
        clade_a.branch_length = max(new_height - heights[name_a], 0.00001)
        clade_b.branch_length = max(new_height - heights[name_b], 0.00001)
        new_clade = Clade(branch_length=0.0, clades=[clade_a, clade_b])

        clusters[new_name] = clusters[name_a] + clusters[name_b]
        clades[new_name] = new_clade
        heights[new_name] = new_height

        del clusters[name_a]
        del clusters[name_b]
        del clades[name_a]
        del clades[name_b]
        del heights[name_a]
        del heights[name_b]

    root_name = next(iter(clades))
    tree = Tree(root=clades[root_name])
    newick = tree.format("newick").strip()

    return newick, tree_to_json(tree.root), distance_matrix, alignment


def tree_to_json(clade: Clade) -> dict:
    return {
        "name": clade.name or "",
        "branch_length": round(float(clade.branch_length or 0.0), 5),
        "children": [tree_to_json(child) for child in clade.clades],
    }


def newick_to_ascii(newick: str) -> str:
    from io import StringIO

    tree = Phylo.read(StringIO(newick), "newick")
    output = StringIO()
    Phylo.draw_ascii(tree, file=output)
    return output.getvalue()
