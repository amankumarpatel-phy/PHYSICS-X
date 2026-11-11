from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from physicsx.phonon import (
    ATOMIC_MASS_UNIT_KG,
    angular_frequency_to_thz,
    diatomic_chain_dispersion,
    phonon_density_of_states,
)

st.set_page_config(page_title="Phonon Lab | PHYSICS-X", page_icon="〰️", layout="wide")
st.title("〰️ Phonon Lab")
st.write(
    "Calculate the acoustic and optical phonon branches of a one-dimensional "
    "diatomic chain using the harmonic nearest-neighbour approximation."
)

with st.sidebar:
    st.header("Lattice parameters")
    lattice_angstrom = st.number_input(
        "Primitive-cell length a (Å)", min_value=0.1, max_value=30.0, value=5.64, step=0.01
    )
    spring_constant = st.number_input(
        "Nearest-neighbour force constant K (N/m)",
        min_value=0.01,
        max_value=2000.0,
        value=20.0,
        step=1.0,
    )
    st.subheader("Atomic masses")
    mass_a_amu = st.number_input(
        "Mass of atom A (u)", min_value=0.1, max_value=300.0, value=22.99, step=0.1
    )
    mass_b_amu = st.number_input(
        "Mass of atom B (u)", min_value=0.1, max_value=300.0, value=35.45, step=0.1
    )
    q_points = st.slider(
        "Wavevector samples", min_value=101, max_value=2001, value=801, step=100
    )

try:
    result = diatomic_chain_dispersion(
        lattice_constant_m=lattice_angstrom * 1e-10,
        force_constant_n_m=spring_constant,
        mass1_kg=mass_a_amu * ATOMIC_MASS_UNIT_KG,
        mass2_kg=mass_b_amu * ATOMIC_MASS_UNIT_KG,
        q_points=q_points,
    )
    acoustic_thz = angular_frequency_to_thz(result["acoustic_omega_rad_s"])
    optical_thz = angular_frequency_to_thz(result["optical_omega_rad_s"])
    dos_frequency, dos_values, dos_edges = phonon_density_of_states(
        result["acoustic_omega_rad_s"],
        result["optical_omega_rad_s"],
        bins=100,
    )
except ValueError as exc:
    st.error(str(exc))
    st.stop()

gamma_index = int(np.argmin(np.abs(result["q_reduced"])))
acoustic_max = float(np.max(acoustic_thz))
optical_gamma = float(optical_thz[gamma_index])
frequency_gap = float(np.min(optical_thz) - np.max(acoustic_thz))

m1, m2, m3 = st.columns(3)
m1.metric("Acoustic branch maximum", f"{acoustic_max:.3f} THz")
m2.metric("Optical frequency at Gamma", f"{optical_gamma:.3f} THz")
m3.metric("Approximate branch gap", f"{frequency_gap:.3f} THz")

dispersion_fig = go.Figure()
dispersion_fig.add_trace(
    go.Scatter(
        x=result["q_reduced"],
        y=acoustic_thz,
        mode="lines",
        name="Acoustic branch",
        line=dict(width=3),
    )
)
dispersion_fig.add_trace(
    go.Scatter(
        x=result["q_reduced"],
        y=optical_thz,
        mode="lines",
        name="Optical branch",
        line=dict(width=3),
    )
)
dispersion_fig.add_vline(x=-1.0, line_dash="dash", annotation_text="-pi/a")
dispersion_fig.add_vline(x=0.0, line_dash="dot", annotation_text="Gamma")
dispersion_fig.add_vline(x=1.0, line_dash="dash", annotation_text="pi/a")
dispersion_fig.update_layout(
    title="Phonon dispersion relation",
    xaxis_title="Reduced wavevector qa/pi",
    yaxis_title="Frequency (THz)",
    template="plotly_dark",
    height=500,
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
)
st.plotly_chart(dispersion_fig, use_container_width=True)

dos_fig = go.Figure(
    go.Bar(
        x=dos_frequency,
        y=dos_values,
        width=np.diff(dos_edges),
        name="Phonon DOS",
    )
)
dos_fig.update_layout(
    title="Estimated phonon density of states",
    xaxis_title="Frequency (THz)",
    yaxis_title="DOS (states / THz / primitive cell)",
    template="plotly_dark",
    height=420,
)
st.plotly_chart(dos_fig, use_container_width=True)

data = pd.DataFrame(
    {
        "q_rad_per_m": result["q_rad_m"],
        "reduced_q_qa_over_pi": result["q_reduced"],
        "acoustic_frequency_THz": acoustic_thz,
        "optical_frequency_THz": optical_thz,
    }
)
st.download_button(
    "Download dispersion data (CSV)",
    data=data.to_csv(index=False).encode("utf-8"),
    file_name="physicsx_diatomic_chain_dispersion.csv",
    mime="text/csv",
)

with st.expander("Model, equations and limitations"):
    st.latex(
        r"\omega_{\pm}^{2}(q)="
        r"K\left(\frac{1}{m_A}+\frac{1}{m_B}\right)"
        r"\pm\sqrt{K^2\left(\frac{1}{m_A}+\frac{1}{m_B}\right)^2"
        r"-\frac{4K^2}{m_A m_B}\sin^2\!\left(\frac{qa}{2}\right)}"
    )
    st.markdown(
        """
        **Model assumptions**

        - One-dimensional, infinite, periodic chain with two atoms in each primitive cell.
        - Alternating atomic masses and identical nearest-neighbour harmonic springs.
        - Small displacements about equilibrium; no anharmonicity or temperature-dependent renormalization.
        - The model includes only longitudinal displacements and short-range spring forces.

        The density of states is a histogram estimate from the sampled wavevectors.
        Its integral is normalized to two vibrational states per primitive cell.

        **Important:** this is not a full three-dimensional lattice-dynamics calculation.
        In particular, it omits transverse branches, long-range Coulomb interactions,
        LO-TO splitting, realistic interatomic force constants and crystal-specific
        Brillouin-zone geometry. For ionic materials, the entered atomic masses alone
        do not make this a quantitative model of a real compound.
        """
    )

st.caption("PHYSICS-X Phonon Lab · Harmonic diatomic-chain model")
