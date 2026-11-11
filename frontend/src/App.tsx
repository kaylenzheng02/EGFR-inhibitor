import { useState } from "react";

function App() {
  const [page, setPage] = useState("Dashboard");

  return (
    <div style={{
      minHeight: "100vh",
      padding: "40px",
      background: "#101820",
      color: "#ffffff",
      fontFamily: "Arial, sans-serif"
    }}>
      <h1>EGFRank</h1>
      <p>Machine Learning-Assisted EGFR Inhibitor Prioritization</p>

      <nav style={{ display: "flex", gap: "12px" }}>
        {["Dashboard", "Compounds", "Model Insights"].map(item => (
          <button
            key={item}
            onClick={() => setPage(item)}
          >
            {item}
          </button>
        ))}
      </nav>

      <main style={{ marginTop: "40px" }}>
        <h2>{page}</h2>
        <p>
          {page === "Dashboard"
            ? "Explore experimental EGFR bioactivity data."
            : page === "Compounds"
            ? "Compare chemical compounds and their activity."
            : "Understand model performance and limitations."}
        </p>
      </main>
    </div>
  );
}

export default App;
