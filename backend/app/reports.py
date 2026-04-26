from __future__ import annotations

import csv
from io import BytesIO, StringIO

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas


def build_csv_summary(analysis: dict, blast_results: list[dict]) -> str:
    buffer = StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["Field", "Value"])
    writer.writerow(["Sequence length", analysis["length"]])
    writer.writerow(["GC percent", analysis["gc_percent"]])
    writer.writerow(["AT percent", analysis["at_percent"]])
    writer.writerow(["Molecular weight Da", analysis["molecular_weight_da"]])
    writer.writerow(["Removed invalid characters", analysis["removed_characters"]])
    writer.writerow(["Longest ORF length", analysis["longest_orf"]["length"] if analysis["longest_orf"] else "None"])
    writer.writerow([])
    writer.writerow(["Organism", "Accession", "Percent identity", "Query coverage", "E-value", "Score"])
    for hit in blast_results:
        writer.writerow(
            [
                hit["organism_name"],
                hit["accession_id"],
                hit["percent_identity"],
                hit["query_coverage"],
                hit["e_value"],
                hit["alignment_score"],
            ]
        )
    return buffer.getvalue()


def build_protein_fasta(analysis: dict) -> str:
    protein = analysis.get("protein_translation", "")
    lines = [">translated_protein_frame_1"]
    lines.extend(protein[index:index + 70] for index in range(0, len(protein), 70))
    return "\n".join(lines) + "\n"


def build_pdf_report(analysis: dict, blast_results: list[dict], newick: str) -> bytes:
    buffer = BytesIO()
    pdf = canvas.Canvas(buffer, pagesize=letter)
    width, height = letter
    y = height - 50

    def line(text: str, step: int = 16) -> None:
        nonlocal y
        if y < 60:
            pdf.showPage()
            y = height - 50
        pdf.drawString(45, y, text[:110])
        y -= step

    pdf.setTitle("DNA Sequence Analysis Report")
    pdf.setFont("Helvetica-Bold", 16)
    line("DNA Sequence Analyzer & Phylogeny Tool", 24)
    pdf.setFont("Helvetica", 10)
    line("Prediction-based educational report. Confirm important results with validated databases.", 22)

    pdf.setFont("Helvetica-Bold", 12)
    line("Basic sequence information")
    pdf.setFont("Helvetica", 10)
    line(f"Length: {analysis['length']} bp")
    line(f"GC: {analysis['gc_percent']}%")
    line(f"AT: {analysis['at_percent']}%")
    line(f"Molecular weight estimate: {analysis['molecular_weight_da']} Da")
    line(f"Invalid/non-DNA characters removed: {analysis['removed_characters']}")

    pdf.setFont("Helvetica-Bold", 12)
    line("Longest ORF")
    pdf.setFont("Helvetica", 10)
    if analysis["longest_orf"]:
        orf = analysis["longest_orf"]
        line(f"Frame {orf['frame']} ({orf['strand']} strand), {orf['start']}-{orf['end']}, {orf['length']} bp")
    else:
        line("No complete ATG-to-stop ORF found.")

    pdf.setFont("Helvetica-Bold", 12)
    line("Top similarity hits")
    pdf.setFont("Helvetica", 10)
    for hit in blast_results[:5]:
        line(
            f"{hit['organism_name']} | {hit['accession_id']} | "
            f"{hit['percent_identity']}% identity | {hit['query_coverage']}% coverage"
        )

    pdf.setFont("Helvetica-Bold", 12)
    line("Phylogenetic tree Newick")
    pdf.setFont("Helvetica", 8)
    for index in range(0, len(newick), 100):
        line(newick[index:index + 100], 12)

    pdf.save()
    return buffer.getvalue()

