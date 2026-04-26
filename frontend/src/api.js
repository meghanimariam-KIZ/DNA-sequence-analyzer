const API_BASE = "http://localhost:8000/api";

export async function analyzeSequence(sequence) {
  const response = await fetch(`${API_BASE}/analyze`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ sequence }),
  });

  if (!response.ok) {
    const error = await response.json().catch(() => ({}));
    throw new Error(error.detail || "Analysis failed.");
  }

  return response.json();
}

export async function analyzeUpload(file) {
  const formData = new FormData();
  formData.append("file", file);

  const response = await fetch(`${API_BASE}/analyze-upload`, {
    method: "POST",
    body: formData,
  });

  if (!response.ok) {
    const error = await response.json().catch(() => ({}));
    throw new Error(error.detail || "Upload analysis failed.");
  }

  return response.json();
}

export async function downloadResult(kind, result) {
  const formData = new FormData();
  formData.append("payload", JSON.stringify(result));

  const response = await fetch(`${API_BASE}/download/${kind}`, {
    method: "POST",
    body: formData,
  });

  if (!response.ok) {
    throw new Error("Download failed.");
  }

  const blob = await response.blob();
  const disposition = response.headers.get("Content-Disposition") || "";
  const filename = disposition.match(/filename=([^;]+)/)?.[1] || `dna-analysis.${kind}`;
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

