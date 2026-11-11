from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from physicsx.quantum import energy_levels_ev, stationary_wavefunction

st.set_page_config(page_title="Quantum Lab | PHYSICS-X", page_icon="⚛️", layout="wide")
st.title("⚛️ Quantum Lab")
st.write("One-dimensional non-relativistic particle in an infinite square well.")

col_a, col_b = st.columns(2)
with col_a:
    width_nm = st.slider("Well width (nm)", min_value=0.1, max_value=5.0, value=1.0, step=0.1)
with col_b:
    particle = st.selectbox("Particle", ["Electron", "Proton"])
mass = {"Electron": 9.1093837139e-31, "Proton": 1.67262192595e-27}[particle]
n_max = st.slider("Number of energy levels", min_value=3, max_value=12, value=6)

width_m = width_nm * 1e-9
levels = energy_levels_ev(width_m, mass, n_max)
x = np.linspace(0.0, width_m, 800)

st.subheader("Energy spectrum")
energy_fig = go.Figure()
for n, energy in enumerate(levels, start=1):
    energy_fig.add_trace(
        go.Scatter(
            x=[0, width_nm],
            y=[energy, energy],
            mode="lines",
            name=f"n = {n}: {energy:.4g} eV",
            line=dict(width=3),
        )
    )
energy_fig.update_layout(
    xaxis_title="Position across well (nm)",
    yaxis_title="Energy (eV)",
    template="plotly_dark",
    height=430,
)
st.plotly_chart(energy_fig, use_container_width=True)

st.subheader("Stationary states")
selected_n = st.slider("Quantum number n", min_value=1, max_value=n_max, value=1)
psi = stationary_wavefunction(width_m, selected_n, x)
probability = np.abs(psi) ** 2
wave_fig = go.Figure()
wave_fig.add_trace(go.Scatter(x=x * 1e9, y=psi, name="ψₙ(x)", line=dict(width=2)))
wave_fig.add_trace(go.Scatter(x=x * 1e9, y=probability / probability.max() * max(np.max(np.abs(psi)), 1e-12), name="|ψₙ(x)|² (scaled)", line=dict(dash="dash")))
wave_fig.update_layout(
    xaxis_title="Position x (nm)",
    yaxis_title="Wavefunction amplitude / scaled probability density",
    template="plotly_dark",
    height=450,
)
st.plotly_chart(wave_fig, use_container_width=True)

selected_energy = levels[selected_n - 1]
m1, m2 = st.columns(2)
m1.metric(f"Energy E{selected_n}", f"{selected_energy:.6g} eV")
m2.metric("Ground-state energy", f"{levels[0]:.6g} eV")

with st.expander("Theory and assumptions"):
    st.latex(r"E_n = \frac{n^2\pi^2\hbar^2}{2mL^2}, \qquad n=1,2,3,\ldots")
    st.latex(r"\psi_n(x)=\sqrt{\frac{2}{L}}\sin\!\left(\frac{n\pi x}{L}\right)")
    st.markdown(
        """
        The potential is zero inside 0 < x < L and infinite outside. The model is
        one-dimensional, non-relativistic and time-independent. Spin, finite barriers,
        interactions, external fields and relativistic corrections are not included.
        The displayed probability curve is scaled only for visual comparison with ψ;
        its physical probability density is |ψ|².
        """
    )
