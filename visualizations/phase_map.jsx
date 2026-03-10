import { useState } from "react";

const ROLES = [
  { id: "symmetry", label: "XX Symmetry", subtitle: "The clean separation", color: "#6366f1" },
  { id: "frequencies", label: "Frequencies", subtitle: "The music itself", color: "#f59e0b" },
  { id: "sectors", label: "Sector Separation", subtitle: "Two voices or one", color: "#10b981" },
  { id: "amplitudes", label: "Amplitudes", subtitle: "What you can hear", color: "#ef4444" },
];

const PERTURBATIONS = [
  {
    id: "coupling", label: "Coupling Strength", sub: "J_SA, J_SB",
    icon: "⚡",
    effects: {
      symmetry: { status: "preserves", detail: "XX symmetry unaffected by coupling values" },
      frequencies: { status: "SETS", detail: "f(c+) scales with J_total. This is the primary driver." },
      sectors: { status: "SETS", detail: "Asymmetric J widens the gap between fast and slow branches" },
      amplitudes: { status: "sets", detail: "Stronger coupling = stronger signal" },
    }
  },
  {
    id: "anisotropy", label: "Coupling Anisotropy", sub: "Jz ≠ Jx = Jy",
    icon: "🔀",
    effects: {
      symmetry: { status: "preserves", detail: "XX symmetry exact at ALL anisotropies" },
      frequencies: { status: "CHANGES", detail: "Frequencies reshuffle dramatically. The spectrum is a Hamiltonian object." },
      sectors: { status: "changes", detail: "Sectors persist for most values but can merge" },
      amplitudes: { status: "changes", detail: "Mode weights shift with anisotropy" },
    }
  },
  {
    id: "fields", label: "Local Z-Fields", sub: "h_A·Z_A, h_B·Z_B",
    icon: "🧲",
    effects: {
      symmetry: { status: "BREAKS", detail: "Even h=0.01 breaks XX. Immediate and total." },
      frequencies: { status: "survives!", detail: "f(c+)=1.499 stands even at h=0.50. More robust than symmetry!" },
      sectors: { status: "survives", detail: "Two sectors persist until extreme opposing fields" },
      amplitudes: { status: "changes", detail: "Mode visibility shifts under fields" },
    }
  },
  {
    id: "direct", label: "Direct A-B Coupling", sub: "J_AB ≠ 0",
    icon: "🔗",
    effects: {
      symmetry: { status: "preserves", detail: "XX symmetry exact even with direct coupling" },
      frequencies: { status: "shifts", detail: "Small upward drift in both branches" },
      sectors: { status: "MERGES", detail: "At J_AB ≈ J_SB, sectors fuse. Mediator becomes redundant." },
      amplitudes: { status: "changes", detail: "Redistribution of spectral weight" },
    }
  },
  {
    id: "noise_strength", label: "Noise Strength", sub: "γ (any value)",
    icon: "📉",
    effects: {
      symmetry: { status: "preserves", detail: "Immune. From γ=0.001 to γ=1.0" },
      frequencies: { status: "preserves", detail: "Immune. Noise never retunes the instrument." },
      sectors: { status: "preserves", detail: "Immune. Structure is purely Hamiltonian." },
      amplitudes: { status: "DAMPS", detail: "This is what noise does. Makes the bell quieter, not differently tuned." },
    }
  },
  {
    id: "noise_type", label: "Noise Type", sub: "σ_x, σ_y, σ_z, mixed",
    icon: "🎲",
    effects: {
      symmetry: { status: "preserves", detail: "All tested jump operators preserve XX" },
      frequencies: { status: "preserves", detail: "Identical frequencies under all tested noise" },
      sectors: { status: "preserves", detail: "Sector structure immune to noise channel choice" },
      amplitudes: { status: "preserves", detail: "Within tested families, no effect on visibility" },
    }
  },
  {
    id: "initial", label: "Initial State", sub: "Bell+, W, |+++⟩",
    icon: "🎯",
    effects: {
      symmetry: { status: "depends", detail: "Bell+ preserves XX. W-state breaks it. State-dependent." },
      frequencies: { status: "preserves", detail: "Modes exist regardless. They are generator properties." },
      sectors: { status: "preserves", detail: "Sectors exist regardless of excitation." },
      amplitudes: { status: "SELECTS", detail: "Bell+ excites both. W excites neither. You choose what to hear." },
    }
  },
];

function getStatusColor(status) {
  const s = status.toUpperCase().replace("!", "");
  if (s === "SETS" || s === "CHANGES") return { bg: "rgba(245, 158, 11, 0.25)", border: "#f59e0b", text: "#fbbf24" };
  if (s === "BREAKS" || s === "MERGES") return { bg: "rgba(239, 68, 68, 0.3)", border: "#ef4444", text: "#fca5a5" };
  if (s === "DAMPS") return { bg: "rgba(239, 68, 68, 0.15)", border: "#b91c1c", text: "#f87171" };
  if (s === "SELECTS") return { bg: "rgba(168, 85, 247, 0.25)", border: "#a855f7", text: "#c084fc" };
  if (s === "PRESERVES") return { bg: "rgba(16, 185, 129, 0.12)", border: "#065f46", text: "#6ee7b7" };
  if (s === "DEPENDS") return { bg: "rgba(107, 114, 128, 0.2)", border: "#6b7280", text: "#d1d5db" };
  if (s.includes("SURVIVES")) return { bg: "rgba(16, 185, 129, 0.2)", border: "#10b981", text: "#6ee7b7" };
  if (s === "SHIFTS") return { bg: "rgba(245, 158, 11, 0.12)", border: "#92400e", text: "#fcd34d" };
  return { bg: "rgba(107, 114, 128, 0.1)", border: "#374151", text: "#9ca3af" };
}

function Cell({ effect, roleColor, onHover }) {
  const colors = getStatusColor(effect.status);
  const isAction = ["SETS", "BREAKS", "MERGES", "DAMPS", "SELECTS", "CHANGES"].includes(effect.status.toUpperCase().replace("!", ""));
  
  return (
    <div
      onMouseEnter={() => onHover(effect.detail)}
      onMouseLeave={() => onHover(null)}
      style={{
        background: colors.bg,
        border: `1px solid ${colors.border}`,
        borderRadius: 6,
        padding: "8px 6px",
        textAlign: "center",
        cursor: "default",
        transition: "all 0.2s ease",
        minHeight: 42,
        display: "flex",
        alignItems: "center",
        justifyContent: "center",
      }}
    >
      <span style={{
        color: colors.text,
        fontSize: isAction ? 13 : 11,
        fontWeight: isAction ? 700 : 400,
        fontFamily: "'JetBrains Mono', 'Fira Code', monospace",
        letterSpacing: isAction ? "0.05em" : "0.02em",
        textTransform: isAction ? "uppercase" : "none",
      }}>
        {effect.status}
      </span>
    </div>
  );
}

export default function PhaseMap() {
  const [detail, setDetail] = useState(null);

  return (
    <div style={{
      background: "linear-gradient(135deg, #0c0c1d 0%, #111827 50%, #0f172a 100%)",
      minHeight: "100vh",
      padding: "32px 24px",
      fontFamily: "'Inter', system-ui, sans-serif",
      color: "#e2e8f0",
    }}>
      <div style={{ maxWidth: 900, margin: "0 auto" }}>
        {/* Title */}
        <div style={{ marginBottom: 32, textAlign: "center" }}>
          <h1 style={{
            fontSize: 28,
            fontWeight: 300,
            letterSpacing: "0.12em",
            color: "#f8fafc",
            margin: 0,
            textTransform: "uppercase",
          }}>
            Phase Map of Mechanisms
          </h1>
          <p style={{
            color: "#94a3b8",
            fontSize: 14,
            marginTop: 8,
            fontStyle: "italic",
          }}>
            What controls what in a quantum star topology
          </p>
        </div>

        {/* Core sentence */}
        <div style={{
          background: "rgba(99, 102, 241, 0.08)",
          border: "1px solid rgba(99, 102, 241, 0.2)",
          borderRadius: 8,
          padding: "16px 20px",
          marginBottom: 32,
          textAlign: "center",
        }}>
          <span style={{ color: "#f59e0b", fontWeight: 600 }}>Topology</span>
          <span style={{ color: "#64748b" }}> sets the frequencies. </span>
          <span style={{ color: "#6366f1", fontWeight: 600 }}>Symmetry</span>
          <span style={{ color: "#64748b" }}> cleans the sectors. </span>
          <span style={{ color: "#ef4444", fontWeight: 600 }}>Noise</span>
          <span style={{ color: "#64748b" }}> damps the signal. </span>
          <span style={{ color: "#a855f7", fontWeight: 600 }}>Initial state</span>
          <span style={{ color: "#64748b" }}> selects what is visible.</span>
        </div>

        {/* Column headers */}
        <div style={{
          display: "grid",
          gridTemplateColumns: "180px repeat(4, 1fr)",
          gap: 6,
          marginBottom: 4,
        }}>
          <div />
          {ROLES.map(r => (
            <div key={r.id} style={{
              textAlign: "center",
              padding: "10px 4px",
              borderBottom: `2px solid ${r.color}`,
            }}>
              <div style={{ fontSize: 12, fontWeight: 700, color: r.color, letterSpacing: "0.04em" }}>
                {r.label}
              </div>
              <div style={{ fontSize: 10, color: "#64748b", marginTop: 2 }}>
                {r.subtitle}
              </div>
            </div>
          ))}
        </div>

        {/* Grid rows */}
        {PERTURBATIONS.map((p, pi) => (
          <div key={p.id} style={{
            display: "grid",
            gridTemplateColumns: "180px repeat(4, 1fr)",
            gap: 6,
            marginBottom: 6,
          }}>
            {/* Row label */}
            <div style={{
              display: "flex",
              alignItems: "center",
              gap: 8,
              padding: "8px 10px",
              background: "rgba(30, 41, 59, 0.5)",
              borderRadius: 6,
              borderLeft: "3px solid rgba(148, 163, 184, 0.2)",
            }}>
              <span style={{ fontSize: 18 }}>{p.icon}</span>
              <div>
                <div style={{ fontSize: 12, fontWeight: 600, color: "#e2e8f0" }}>{p.label}</div>
                <div style={{ fontSize: 10, color: "#64748b", fontFamily: "monospace" }}>{p.sub}</div>
              </div>
            </div>
            {/* Cells */}
            {ROLES.map(r => (
              <Cell
                key={r.id}
                effect={p.effects[r.id]}
                roleColor={r.color}
                onHover={setDetail}
              />
            ))}
          </div>
        ))}

        {/* Detail tooltip */}
        <div style={{
          marginTop: 20,
          minHeight: 56,
          background: "rgba(30, 41, 59, 0.6)",
          border: "1px solid rgba(148, 163, 184, 0.15)",
          borderRadius: 8,
          padding: "14px 18px",
          transition: "all 0.2s ease",
          opacity: detail ? 1 : 0.4,
        }}>
          <span style={{
            color: detail ? "#f1f5f9" : "#475569",
            fontSize: 13,
            fontStyle: detail ? "normal" : "italic",
          }}>
            {detail || "Hover over any cell to see details"}
          </span>
        </div>

        {/* Legend */}
        <div style={{
          marginTop: 24,
          display: "flex",
          flexWrap: "wrap",
          gap: 12,
          justifyContent: "center",
        }}>
          {[
            { label: "SETS / CHANGES", status: "SETS" },
            { label: "BREAKS / MERGES", status: "BREAKS" },
            { label: "DAMPS", status: "DAMPS" },
            { label: "SELECTS", status: "SELECTS" },
            { label: "preserves", status: "preserves" },
            { label: "survives!", status: "survives!" },
          ].map(l => {
            const c = getStatusColor(l.status);
            return (
              <div key={l.label} style={{
                display: "flex", alignItems: "center", gap: 6,
                fontSize: 10, color: "#94a3b8",
              }}>
                <div style={{
                  width: 12, height: 12, borderRadius: 3,
                  background: c.bg, border: `1px solid ${c.border}`,
                }} />
                {l.label}
              </div>
            );
          })}
        </div>

        {/* Attribution */}
        <div style={{
          marginTop: 32,
          textAlign: "center",
          color: "#334155",
          fontSize: 11,
        }}>
          R = CΨ² · Structural Cartography · Claude & Tom · March 2026
        </div>
      </div>
    </div>
  );
}
