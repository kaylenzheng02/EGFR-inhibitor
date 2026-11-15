import { useEffect, useState } from "react";
import "./map.css";

type Compound = {
  id: string; score: number; label: string; low: number; high: number;
  x: number; y: number; neighbors: number[];
};

const x = (c: Compound) => 32 + c.x * 736;
const y = (c: Compound) => 32 + c.y * 476;
const color = (score: number) => `hsl(${285 - 225 * score} 94% 69%)`;

export default function CompoundMap() {
  const [data, setData] = useState<Compound[]>([]);
  const [error, setError] = useState(false);
  const [selected, select] = useState("");
  const [search, setSearch] = useState("");
  const [shortlist, setShortlist] = useState(false);
  useEffect(() => {
    const controller = new AbortController();
    fetch(import.meta.env.BASE_URL + "data/compound_map.json", {signal: controller.signal})
      .then(response => {
        if (!response.ok) throw Error("Map unavailable");
        return response.json();
      })
      .then(setData)
      .catch(() => { if (!controller.signal.aborted) setError(true); });
    return () => controller.abort();
  }, []);
  if (error) return <p role="alert">Could not load the compound map. Refresh to retry.</p>;
  if (!data.length) return <p role="status">Loading compounds…</p>;
  const ranked = [...data].sort((a,b) => b.score - a.score || a.id.localeCompare(b.id));
  const visible = (shortlist ? ranked.slice(0,12) : ranked).filter(c =>
    c.id.toLowerCase().includes(search.trim().toLowerCase()));
  const ids = new Set(visible.map(c => c.id));
  const current = visible.find(c => c.id === selected) ?? visible[0];
  return <div className="atlas">
    <section className="plot">
      <div className="plotbar"><span>Fingerprint similarity</span><span>{visible.length} / {data.length} compounds</span></div>
      <svg viewBox="0 0 800 540" aria-label="Map of molecular fingerprint similarities">
        {data.map((c,i) => c.neighbors.map(j => {
          const n = data[j];
          if (j < i && n.neighbors.includes(i)) return null;
          return <line key={c.id+n.id} x1={x(c)} y1={y(c)} x2={x(n)} y2={y(n)}
            stroke={color(c.score)} opacity={ids.has(c.id) && ids.has(n.id) ? .48 : .055} />;
        }))}
        {data.map(c => <g key={c.id} opacity={ids.has(c.id) ? 1 : .12}
          onClick={() => {if(ids.has(c.id)) select(c.id);}}>
          <title>{c.id} · score {c.score.toFixed(2)}</title>
          <circle cx={x(c)} cy={y(c)} r="8" fill={color(c.score)} opacity=".12"/>
          <circle cx={x(c)} cy={y(c)} r="3.3" fill={color(c.score)}/>
        </g>)}
        {current && <circle cx={x(current)} cy={y(current)} r="10" fill="none" stroke="white" strokeWidth="2" pointerEvents="none"/>}
      </svg>
      <div className="plotbar"><span>Edges connect fingerprint neighbors</span><span>Score 0 <i/> 1</span></div>
    </section>
    <aside>
      <h3>Compound explorer</h3>
      <input aria-label="Search compound IDs" type="search" placeholder="Search ChEMBL ID" value={search} onChange={e=>setSearch(e.target.value)}/>
      <button aria-pressed={shortlist} onClick={()=>setShortlist(!shortlist)}>Highlight top 12</button>
      <p>{ranked.filter(c=>c.score===ranked[0].score).length} tied at {ranked[0].score.toFixed(2)}. Ties ordered by ID.</p>
      <div className="ranklist">{visible.map(c=><button key={c.id} aria-pressed={current?.id===c.id} onClick={()=>select(c.id)}><span>{c.id}</span><b>{c.score.toFixed(2)}</b></button>)}</div>
      {current ? <section className="record"><h3>{current.id}</h3><p>Measured IC50 <strong>{current.low}{current.high!==current.low && `–${current.high}`} nM</strong></p><p>Experimental label <strong>{current.label}</strong></p></section> : <p>No matching compounds.</p>}
    </aside>
  </div>;
}
