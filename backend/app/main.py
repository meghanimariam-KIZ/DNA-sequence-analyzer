from __future__ import annotations

import json
from typing import Literal

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response
from pydantic import BaseModel

from .bio_utils import analyze_sequence
from .blast import runBlastSearch
from .phylo import build_upgma_tree, newick_to_ascii
from .reports import build_csv_summary, build_pdf_report, build_protein_fasta


app = FastAPI(
    title="DNA Sequence Analyzer & Phylogeny Tool",
    description="Educational DNA analysis, demo similarity search, and phylogeny API.",
    version="1.0.0",
)
@app.get("/")
def home():
    return {"message": "DNA Analyzer Backend is LIVE"}
    
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://127.0.0.1:5173",
        "https://dna-sequence-analyzer-five.vercel.app/",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class AnalyzeRequest(BaseModel):
    sequence: str


def _full_analysis(sequence_text: str) -> dict:
    analysis = analyze_sequence(sequence_text)
    blast_results = runBlastSearch(analysis["cleaned_sequence"])

    phylo_sequences = {"Query": analysis["cleaned_sequence"]}
    for hit in blast_results[:5]:
        label = f"{hit['accession_id']} {hit['organism_name']}"
        phylo_sequences[label] = hit["sequence"]

    newick, tree_json, distance_matrix, alignment = build_upgma_tree(phylo_sequences)
    for hit in blast_results:
        hit.pop("sequence", None)

    return {
        "analysis": analysis,
        "similarity_results": blast_results,
        "phylogeny": {
            "method": "Simple star MSA followed by UPGMA with p-distance",
            "newick": newick,
            "ascii_tree": newick_to_ascii(newick),
            "tree": tree_json,
            "distance_matrix": distance_matrix,
            "alignment_preview": {
                name: aligned_sequence[:180] for name, aligned_sequence in alignment.items()
            },
        },
        "scientific_note": (
            "This result is prediction-based. Demo similarity search uses a small local "
            "reference database and should not be interpreted as verified organism identification."
        ),
    }


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/api/analyze")
def analyze(request: AnalyzeRequest) -> dict:
    try:
        return _full_analysis(request.sequence)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/analyze-upload")
async def analyze_upload(file: UploadFile = File(...)) -> dict:
    content = await file.read()
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HTTPException(status_code=400, detail="Uploaded file must be plain text or FASTA.") from exc

    try:
        return _full_analysis(text)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/download/{kind}")
def download(
    kind: Literal["pdf", "csv", "newick", "protein"],
    payload: str = Form(...),
) -> Response:
    try:
        data = json.loads(payload)
        analysis = data["analysis"]
        blast_results = data["similarity_results"]
        newick = data["phylogeny"]["newick"]
    except (KeyError, json.JSONDecodeError) as exc:
        raise HTTPException(status_code=400, detail="Download payload is incomplete.") from exc

    if kind == "csv":
        content = build_csv_summary(analysis, blast_results)
        return Response(
            content,
            media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=sequence_summary.csv"},
        )

    if kind == "newick":
        return Response(
            newick + "\n",
            media_type="text/plain",
            headers={"Content-Disposition": "attachment; filename=phylogenetic_tree.nwk"},
        )

    if kind == "protein":
        content = build_protein_fasta(analysis)
        return Response(
            content,
            media_type="text/plain",
            headers={"Content-Disposition": "attachment; filename=translated_protein.fasta"},
        )

    pdf = build_pdf_report(analysis, blast_results, newick)
    return Response(
        pdf,
        media_type="application/pdf",
        headers={"Content-Disposition": "attachment; filename=dna_analysis_report.pdf"},
    )
