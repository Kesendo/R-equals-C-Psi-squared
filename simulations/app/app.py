"""
Five Regulator Simulator
Interactive 3-qubit open Heisenberg star (S-A-B) with real-time control.
"""
import numpy as np
import streamlit as st
import plotly.graph_objects as go
from physics import run_simulation, INITIAL_STATES

st.set_page_config(
    page_title="Five Regulator Simulator",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded",
)

# === Colors ===
C_PLUS = "#00D4FF"
C_MINUS = "#FF6B35"
C_CONC = "#A855F7"
C_PSI = "#22C55E"
C_CPSI = "#FACC15"

# === Chart layout helper ===
_DARK_MARGIN = dict(l=50, r=20, t=40, b=40)


def dark_layout(fig, title, height=250, **kwargs):
    """Apply consistent dark theme layout to a Plotly figure."""
    kwargs.setdefault('margin', _DARK_MARGIN)
    fig.update_layout(
        title=title,
        template="plotly_dark",
        height=height,
        **kwargs,
    )


# === Sidebar: 5 Regulators ===
with st.sidebar:
    st.title("⚛️ Five Regulators")

    st.header("1 · Topology")
    J_SA = st.slider("J_SA (S–A)", 0.1, 5.0, 1.0, 0.1)
    J_SB = st.slider("J_SB (S–B)", 0.1, 5.0, 2.0, 0.1)

    st.header("2 · Symmetry")
    xy_ratio = st.slider("XY ratio", 0.0, 1.0, 1.0, 0.05,
                          help="1 = isotropic Heisenberg, 0 = pure ZZ (Ising)")

    st.header("3 · Noise strength")
    gamma_exp = st.slider("log₁₀(γ)", -3.0, 0.0, -1.3, 0.05)
    gamma = 10 ** gamma_exp
    st.caption(f"γ = {gamma:.4f}")

    st.header("4 · Initial state")
    state_name = st.selectbox("State", list(INITIAL_STATES.keys()))

    st.header("5 · Bath geometry")
    eta = st.slider("η  (0 = local, 1 = correlated)", 0.0, 1.0, 0.0, 0.05)
    phi = st.slider("φ  (0 = ZZ, π/2 = XX)", 0.0, np.pi / 2, 0.0, 0.05,
                     format="%.2f")
    if phi > 0:
        st.caption(f"φ = {phi:.2f} rad  ({np.degrees(phi):.0f}°)")

    st.divider()
    t_max = st.slider("Simulation time", 5.0, 60.0, 20.0, 5.0)


# === Run simulation (cached) ===
@st.cache_data
def cached_simulation(J_SA, J_SB, xy_ratio, gamma, eta, phi, state_name, t_max):
    return run_simulation(J_SA, J_SB, xy_ratio, gamma, eta, phi, state_name,
                          t_max=t_max)


res = cached_simulation(J_SA, J_SB, xy_ratio, gamma, eta, phi, state_name, t_max)

# === Title ===
st.markdown("## ⚛️ Five Regulator Simulator")
st.caption("3-qubit open Heisenberg star  S(0)–A(1)–B(2)  ·  Lindblad RK4")

# === Main plot: c+(t) and c-(t) ===
fig = go.Figure()
fig.add_trace(go.Scatter(x=res['times'], y=res['c_plus'],
                          mode='lines', name='c₊(t)',
                          line=dict(color=C_PLUS, width=2)))
fig.add_trace(go.Scatter(x=res['times'], y=res['c_minus'],
                          mode='lines', name='c₋(t)',
                          line=dict(color=C_MINUS, width=2)))
dark_layout(fig, "Sector Dynamics", height=400,
            xaxis_title="Time", yaxis_title="Amplitude",
            margin=dict(l=60, r=30, t=50, b=50),
            legend=dict(orientation="h", yanchor="bottom", y=1.02,
                        xanchor="right", x=1))
st.plotly_chart(fig, use_container_width=True)

# === Optional: C*Psi plot ===
with st.expander("C·Ψ dynamics", expanded=False):
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=res['times'], y=res['concurrence'],
                               mode='lines', name='C (concurrence)',
                               line=dict(color=C_CONC, width=1.5)))
    fig2.add_trace(go.Scatter(x=res['times'], y=res['psi'],
                               mode='lines', name='Ψ (coherence)',
                               line=dict(color=C_PSI, width=1.5)))
    fig2.add_trace(go.Scatter(x=res['times'], y=res['cpsi'],
                               mode='lines', name='C·Ψ',
                               line=dict(color=C_CPSI, width=2)))
    fig2.add_hline(y=0.25, line_dash="dot", line_color="white",
                    annotation_text="¼ boundary", annotation_font_color="white")
    dark_layout(fig2, "Concurrence, Coherence & C·Ψ", height=350,
                xaxis_title="Time", yaxis_title="Value",
                margin=dict(l=60, r=30, t=50, b=50))
    st.plotly_chart(fig2, use_container_width=True)

# === Indicators ===
col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("### Frequencies")
    m1, m2 = st.columns(2)
    m1.metric("f(c₊)", f"{res['f_cp']:.3f}")
    m2.metric("f(c₋)", f"{res['f_cm']:.3f}")

    # FFT plot
    f_max = min(5.0, res['fft_freqs'][-1]) if len(res['fft_freqs']) > 0 else 5.0
    mask = res['fft_freqs'] <= f_max
    fig_fft = go.Figure()
    fig_fft.add_trace(go.Scatter(x=res['fft_freqs'][mask], y=res['fft_cp'][mask],
                                  mode='lines', name='c₊',
                                  line=dict(color=C_PLUS, width=1.5)))
    fig_fft.add_trace(go.Scatter(x=res['fft_freqs'][mask], y=res['fft_cm'][mask],
                                  mode='lines', name='c₋',
                                  line=dict(color=C_MINUS, width=1.5)))
    dark_layout(fig_fft, "FFT spectrum",
                xaxis_title="Frequency", yaxis_title="Amplitude",
                showlegend=True, legend=dict(font=dict(size=10)))
    st.plotly_chart(fig_fft, use_container_width=True)

with col2:
    st.markdown("### Amplitudes")
    ratio_str = f"{res['amp_ratio']:.2f}" if res['amp_ratio'] < 100 else "∞"
    st.metric("A₊ / A₋", ratio_str)

    fig_bar = go.Figure()
    fig_bar.add_trace(go.Bar(
        x=["A₊", "A₋"],
        y=[res['a_cp'], res['a_cm']],
        marker_color=[C_PLUS, C_MINUS],
        text=[f"{res['a_cp']:.4f}", f"{res['a_cm']:.4f}"],
        textposition='outside',
    ))
    dark_layout(fig_bar, "Peak amplitudes", yaxis_title="Amplitude")
    st.plotly_chart(fig_bar, use_container_width=True)

with col3:
    st.markdown("### Structure")
    if res['xx_sym']:
        st.success(f"XX symmetry ✓  ({res['xx_final']:.1e})")
    else:
        st.error(f"XX symmetry ✗  ({res['xx_final']:.2e})")

    st.metric("Skeleton fraction", f"{res['skeleton'] * 100:.1f} %")

    # XX commutator over time
    fig_xx = go.Figure()
    fig_xx.add_trace(go.Scatter(x=res['times'], y=res['xx_comms'],
                                 mode='lines', name='‖[ρ,XX]‖',
                                 line=dict(color="#F472B6", width=1.5)))
    fig_xx.add_hline(y=1e-6, line_dash="dot", line_color="gray",
                      annotation_text="threshold")
    dark_layout(fig_xx, "XX commutator", xaxis_title="Time", yaxis_type="log")
    st.plotly_chart(fig_xx, use_container_width=True)

# === Phase map summary ===
st.divider()
dom = "c₊ dominates" if res['a_cp'] > res['a_cm'] else "c₋ dominates"
sym_txt = "preserved" if res['xx_sym'] else f"broken ({res['xx_final']:.1e})"
cpsi_max = np.max(res['cpsi'])
cpsi_status = f"max C·Ψ = {cpsi_max:.3f}" + (" ≥ ¼ !" if cpsi_max >= 0.25 else " < ¼")

st.markdown(
    f"**Topology:** f(c₊) = {res['f_cp']:.3f}, f(c₋) = {res['f_cm']:.3f} &ensp;|&ensp; "
    f"**Noise:** γ = {gamma:.4f} &ensp;|&ensp; "
    f"**Bath:** {dom} ({ratio_str}) &ensp;|&ensp; "
    f"**Symmetry:** {sym_txt} &ensp;|&ensp; "
    f"**Skeleton:** {res['skeleton'] * 100:.1f}% &ensp;|&ensp; "
    f"**{cpsi_status}**"
)
