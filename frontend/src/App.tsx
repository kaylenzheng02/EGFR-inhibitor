import { useEffect, useState } from "react";

type Measurement = {
  activity_id: string;
  molecule_chembl_id: string;
  assay_chembl_id: string;
  standard_value: string;
  pIC50: number;
};

const endpoints: Record<string, string> = {
  CHEMBL674637: "EGFR kinase tyrosine phosphorylation",
  CHEMBL621151: "EGFR autophosphorylation in A431 cells",
  CHEMBL615325: "Ligand-induced proliferation in 3T3 cells",
  CHEMBL674643: "EGF-dependent receptor autophosphorylation",
  CHEMBL677833: "Kinase isolated from A431 cells",
  CHEMBL679674: "EGFR inhibition (readout unspecified)",
};

export default function App() {
  const [page, setPage] = useState("Dashboard");
  const [rows, setRows] = useState<Measurement[] | null>(null);
  const [failed, setFailed] = useState(false);
  const [search, setSearch] = useState("");

  useEffect(() => {
    const controller = new AbortController();
    async function load() {
      try {
        const response = await fetch(import.meta.env.BASE_URL + "data/egfr_sample.json", { signal: controller.signal });
        if (!response.ok) throw new Error("Dataset unavailable");
        setRows(await response.json());
      } catch {
        if (!controller.signal.aborted) setFailed(true);
      }
    }
    void load();
    return () => controller.abort();
  }, []);

  const filtered = (rows ?? []).filter(row =>
    row.molecule_chembl_id.toLowerCase().includes(search.trim().toLowerCase())
  );

  return <div className="shell">
    <header><h1>EGFRank</h1><p>Machine Learning-Assisted EGFR Inhibitor Prioritization</p></header>
    <nav aria-label="Sections">{["Dashboard", "Compounds", "Model Insights"].map(name =>
      <button key={name} aria-pressed={page === name} onClick={() => setPage(name)}>{name}</button>
    )}</nav>
    <main><h2>{page}</h2>
      {page === "Model Insights" ? <p>Model training is planned. The explorer shows experimental measurements.</p>
        : failed ? <p role="alert">Could not load the sample. Run the preprocessing script and refresh this page.</p>
        : rows === null ? <p role="status">Loading measurements…</p>
        : page === "Dashboard" ? <>
          <p>Screened EGFR measurements from ChEMBL</p>
          <dl>{Object.entries({
            Measurements: rows.length,
            Compounds: new Set(rows.map(row => row.molecule_chembl_id)).size,
            Assays: new Set(rows.map(row => row.assay_chembl_id)).size,
          }).map(([label, value]) => <div key={label}><dt>{label}</dt><dd>{value}</dd></div>)}</dl>
        </> : <>
          <label>Search by compound ID<input type="search" placeholder="e.g. CHEMBL68920" value={search} onChange={e => setSearch(e.target.value)} /></label>
          <p aria-live="polite">{filtered.length} matching measurements</p>
          {filtered.length > 0 && <div className="table-wrap"><table>
            <thead><tr>{["Compound ID", "Assay ID", "What the assay measured", "IC50 (nM)", "pIC50"].map(label => <th key={label} scope="col">{label}</th>)}</tr></thead>
            <tbody>{filtered.map(row => <tr key={row.activity_id}>
              <td>{row.molecule_chembl_id}</td><td>{row.assay_chembl_id}</td>
              <td>{endpoints[row.assay_chembl_id] ?? "Assay awaiting review"}</td>
              <td>{row.standard_value}</td><td>{row.pIC50.toFixed(4)}</td>
            </tr>)}</tbody>
          </table></div>}
        </>}
    </main>
  </div>;
}
