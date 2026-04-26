import { useMemo, useState } from "react";
import Plot from "react-plotly.js";
import {
  Activity,
  Download,
  Dna,
  FileText,
  FlaskConical,
  GitBranch,
  Upload,
} from "lucide-react";
import { analyzeSequence, analyzeUpload, downloadResult } from "./api";

const tabs = [
  "Basic Sequence Info",
  "ORF & Translation",
  "Similarity Search",
  "Phylogenetic Tree",
  "Downloads",
];

const exampleSequence = `>example_ecoli_like_sequence
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGA`;

function Metric({ label, value }) {
  return (
    <div className="metric">
      <span>{label}</span>
      <strong>{value}</strong>
    </div>
  );
}

function SequenceBlock({ title, value }) {
  return (
    <section className="sequence-block">
      <h3>{title}</h3>
      <pre>{value || "No sequence available."}</pre>
    </section>
  );
}

function TreePlot({ tree }) {
  const data = useMemo(() => {
    if (!tree) return null;
    const x = [];
    const y = [];
    const text = [];
    const edgesX = [];
    const edgesY = [];
    let leafIndex = 0;

    function visit(node, depth = 0) {
      const children = node.children || [];
      let currentY;
      if (children.length === 0) {
        currentY = leafIndex;
        leafIndex += 1;
      } else {
        const childPoints = children.map((child) => visit(child, depth + (child.branch_length || 0.1)));
        currentY = childPoints.reduce((sum, point) => sum + point.y, 0) / childPoints.length;
        childPoints.forEach((point) => {
          edgesX.push(depth, point.x, null);
          edgesY.push(currentY, point.y, null);
        });
      }
      x.push(depth);
      y.push(currentY);
      text.push(node.name || "internal node");
      return { x: depth, y: currentY };
    }

    visit(tree, 0);

    return {
      nodes: { x, y, text },
      edges: { x: edgesX, y: edgesY },
    };
  }, [tree]);

  if (!data) {
    return <p className="muted">Run an analysis to display a tree.</p>;
  }

  return (
    <Plot
      data={[
        {
          x: data.edges.x,
          y: data.edges.y,
          mode: "lines",
          type: "scatter",
          hoverinfo: "skip",
          line: { color: "#334155", width: 2 },
        },
        {
          x: data.nodes.x,
          y: data.nodes.y,
          text: data.nodes.text,
          mode: "markers+text",
          type: "scatter",
          textposition: "right center",
          marker: { color: "#0f766e", size: 9 },
        },
      ]}
      layout={{
        autosize: true,
        height: 440,
        margin: { l: 20, r: 180, t: 20, b: 40 },
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        xaxis: { title: "Genetic distance", zeroline: false },
        yaxis: { visible: false },
        showlegend: false,
      }}
      useResizeHandler
      className="plot"
      config={{ displayModeBar: false, responsive: true }}
    />
  );
}

export default function App() {
  const [sequence, setSequence] = useState(exampleSequence);
  const [activeTab, setActiveTab] = useState(tabs[0]);
  const [result, setResult] = useState(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [fileName, setFileName] = useState("");

  async function runAnalysis() {
    setLoading(true);
    setError("");
    try {
      const data = await analyzeSequence(sequence);
      setResult(data);
      setActiveTab(tabs[0]);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  }

  async function handleUpload(event) {
    const file = event.target.files?.[0];
    if (!file) return;
    setLoading(true);
    setError("");
    setFileName(file.name);
    try {
      const text = await file.text();
      setSequence(text);
      const data = await analyzeUpload(file);
      setResult(data);
      setActiveTab(tabs[0]);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  }

  const analysis = result?.analysis;

  return (
    <main className="app-shell">
      <header className="topbar">
        <div>
          <h1>DNA Sequence Analyzer & Phylogeny Tool</h1>
          <p>Paste or upload DNA, then explore sequence statistics, ORFs, demo similarity hits, and a predicted tree.</p>
        </div>
        <div className="status-pill">
          <Dna size={18} />
          Educational demo mode
        </div>
      </header>

      <section className="notice">
        <FlaskConical size={20} />
        <span>{result?.scientific_note || "Predictions depend on the reference database used. Confirm important identifications with validated BLAST/database results."}</span>
      </section>

      <section className="workspace">
        <aside className="input-panel">
          <div className="panel-title">
            <FileText size={20} />
            <h2>Input DNA Sequence</h2>
          </div>
          <textarea
            value={sequence}
            onChange={(event) => setSequence(event.target.value)}
            spellCheck="false"
            aria-label="DNA sequence input"
          />
          <div className="actions">
            <label className="upload-button">
              <Upload size={18} />
              Upload FASTA
              <input type="file" accept=".fa,.fasta,.txt" onChange={handleUpload} />
            </label>
            <button className="primary-button" onClick={runAnalysis} disabled={loading}>
              <Activity size={18} />
              {loading ? "Analyzing..." : "Analyze"}
            </button>
          </div>
          {fileName && <p className="muted">Uploaded: {fileName}</p>}
          {error && <p className="error">{error}</p>}
        </aside>

        <section className="results-panel">
          <nav className="tabs" aria-label="Result tabs">
            {tabs.map((tab) => (
              <button
                key={tab}
                className={activeTab === tab ? "active" : ""}
                onClick={() => setActiveTab(tab)}
              >
                {tab}
              </button>
            ))}
          </nav>

          {!result && (
            <div className="empty-state">
              <Dna size={42} />
              <h2>Ready for analysis</h2>
              <p>Use the built-in example, paste your own sequence, or upload a FASTA file.</p>
            </div>
          )}

          {result && activeTab === "Basic Sequence Info" && (
            <div className="tab-content">
              <div className="metrics-grid">
                <Metric label="Sequence length" value={`${analysis.length} bp`} />
                <Metric label="GC content" value={`${analysis.gc_percent}%`} />
                <Metric label="AT content" value={`${analysis.at_percent}%`} />
                <Metric label="Molecular weight" value={`${analysis.molecular_weight_da} Da`} />
                <Metric label="Removed characters" value={analysis.removed_characters} />
                <Metric label="ORFs found" value={analysis.orfs.length} />
              </div>
              <SequenceBlock title="Cleaned DNA Sequence" value={analysis.cleaned_sequence} />
              <SequenceBlock title="Reverse Complement" value={analysis.reverse_complement} />
              <SequenceBlock title="Transcribed RNA" value={analysis.rna_sequence} />
            </div>
          )}

          {result && activeTab === "ORF & Translation" && (
            <div className="tab-content">
              <SequenceBlock title="Protein Translation, Frame +1" value={analysis.protein_translation} />
              <div className="two-column">
                <section>
                  <h3>Longest ORF</h3>
                  {analysis.longest_orf ? (
                    <table>
                      <tbody>
                        <tr><th>Frame</th><td>{analysis.longest_orf.frame}</td></tr>
                        <tr><th>Strand</th><td>{analysis.longest_orf.strand}</td></tr>
                        <tr><th>Coordinates</th><td>{analysis.longest_orf.start}-{analysis.longest_orf.end}</td></tr>
                        <tr><th>Length</th><td>{analysis.longest_orf.length} bp</td></tr>
                        <tr><th>Stop codon</th><td>{analysis.longest_orf.stop_codon}</td></tr>
                      </tbody>
                    </table>
                  ) : (
                    <p className="muted">No complete ATG-to-stop ORF found.</p>
                  )}
                </section>
                <section>
                  <h3>Reading Frames</h3>
                  <table>
                    <thead>
                      <tr><th>Frame</th><th>Strand</th><th>Protein length</th><th>ORFs</th></tr>
                    </thead>
                    <tbody>
                      {analysis.reading_frames.map((frame) => (
                        <tr key={`${frame.strand}-${frame.frame}`}>
                          <td>{frame.frame}</td>
                          <td>{frame.strand}</td>
                          <td>{frame.protein_length}</td>
                          <td>{frame.orf_count}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </section>
              </div>
              <section>
                <h3>Codon Usage</h3>
                <div className="codon-grid">
                  {Object.entries(analysis.codon_usage).map(([codon, count]) => (
                    <span key={codon}><strong>{codon}</strong> {count}</span>
                  ))}
                </div>
              </section>
            </div>
          )}

          {result && activeTab === "Similarity Search" && (
            <div className="tab-content">
              <h3>Closest demo database matches</h3>
              <table>
                <thead>
                  <tr>
                    <th>Organism</th>
                    <th>Accession</th>
                    <th>% Identity</th>
                    <th>Coverage</th>
                    <th>E-value</th>
                    <th>Score</th>
                  </tr>
                </thead>
                <tbody>
                  {result.similarity_results.map((hit) => (
                    <tr key={hit.accession_id}>
                      <td>{hit.organism_name}</td>
                      <td>{hit.accession_id}</td>
                      <td>{hit.percent_identity}</td>
                      <td>{hit.query_coverage}</td>
                      <td>{hit.e_value}</td>
                      <td>{hit.alignment_score}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}

          {result && activeTab === "Phylogenetic Tree" && (
            <div className="tab-content">
              <div className="section-heading">
                <GitBranch size={20} />
                <h3>{result.phylogeny.method}</h3>
              </div>
              <TreePlot tree={result.phylogeny.tree} />
              <section>
                <h3>Multiple Sequence Alignment Preview</h3>
                <div className="alignment-list">
                  {Object.entries(result.phylogeny.alignment_preview).map(([name, aligned]) => (
                    <div key={name}>
                      <strong>{name}</strong>
                      <pre>{aligned}</pre>
                    </div>
                  ))}
                </div>
              </section>
              <SequenceBlock title="Newick Tree" value={result.phylogeny.newick} />
              <SequenceBlock title="ASCII Tree" value={result.phylogeny.ascii_tree} />
            </div>
          )}

          {result && activeTab === "Downloads" && (
            <div className="downloads-grid">
              {[
                ["pdf", "Analysis report PDF"],
                ["csv", "Sequence summary CSV"],
                ["newick", "Phylogenetic tree Newick"],
                ["protein", "Translated protein FASTA"],
              ].map(([kind, label]) => (
                <button key={kind} onClick={() => downloadResult(kind, result)}>
                  <Download size={20} />
                  {label}
                </button>
              ))}
            </div>
          )}
        </section>
      </section>
    </main>
  );
}
