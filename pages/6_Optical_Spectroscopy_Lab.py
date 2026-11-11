from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from physicsx.optics import (
    LorentzOscillator,
    drude_lorentz_dielectric,
    kramers_kronig_eps1_finite_window,
    normal_incidence_slab_rt,
    optical_constants,
    optical_conductivity_s_m,
)

st.set_page_config(
    page_title="Optical Spectroscopy Lab | PHYSICS-X",
    page_icon="🌈",
    layout="wide",
)

st.title("🌈 Optical Spectroscopy Laboratory")
st.write(
    "Explore dielectric response, complex refractive index, absorption, "
    "thin-film R/T/A, optical conductivity and a finite-window Kramers–Kronig check."
)
st.caption("Model: isotropic Drude–Lorentz dielectric function · photon energies in eV")

with st.sidebar:
    st.header("Spectral range")
    energy_min = st.number_input(
        "Minimum photon energy (eV)", min_value=0.01, max_value=20.0, value=0.20, step=0.05
    )
    energy_max = st.number_input(
        "Maximum photon energy (eV)", min_value=0.02, max_value=40.0, value=6.00, step=0.10
    )
    point_count = st.slider("Spectral points", min_value=301, max_value=2401, value=1001, step=100)

    st.header("Background and Drude response")
    eps_inf = st.number_input("High-frequency dielectric constant ε∞", min_value=0.01, max_value=30.0, value=1.0, step=0.1)
    plasma_energy = st.number_input("Plasma energy Ep (eV)", min_value=0.0, max_value=30.0, value=9.0, step=0.1)
    drude_gamma = st.number_input("Drude damping ΓD (eV)", min_value=0.001, max_value=10.0, value=0.35, step=0.01)

    st.header("Bound-electron resonances")
    oscillator_count = st.selectbox("Number of Lorentz oscillators", [0, 1, 2, 3], index=1)
    oscillator_defaults = [
        (2.5, 0.35, 0.80),
        (4.2, 0.55, 0.35),
        (5.4, 0.70, 0.20),
    ]
    oscillators = []
    for i in range(oscillator_count):
        e0_default, gamma_default, strength_default = oscillator_defaults[i]
        with st.expander(f"Lorentz oscillator {i + 1}", expanded=(i == 0)):
            e0 = st.number_input(
                f"Resonance E₀{i + 1} (eV)", min_value=0.01, max_value=30.0,
                value=e0_default, step=0.05, key=f"opt_e0_{i}"
            )
            gamma = st.number_input(
                f"Damping Γ{i + 1} (eV)", min_value=0.001, max_value=10.0,
                value=gamma_default, step=0.01, key=f"opt_gamma_{i}"
            )
            strength = st.number_input(
                f"Oscillator strength f{i + 1}", min_value=0.0, max_value=10.0,
                value=strength_default, step=0.05, key=f"opt_strength_{i}"
            )
            oscillators.append(LorentzOscillator(e0, gamma, strength))

    st.header("Thin-film stack")
    film_thickness_nm = st.number_input(
        "Film thickness (nm)", min_value=0.0, max_value=10000.0, value=250.0, step=10.0
    )
    incident_index = st.number_input(
        "Incident-medium refractive index", min_value=0.1, max_value=5.0, value=1.0, step=0.05
    )
    substrate_index = st.number_input(
        "Substrate refractive index", min_value=0.1, max_value=5.0, value=1.5, step=0.05
    )

if energy_max <= energy_min:
    st.error("Maximum photon energy must be greater than the minimum.")
    st.stop()

energy = np.linspace(energy_min, energy_max, point_count)
try:
    epsilon = drude_lorentz_dielectric(
        energy, eps_inf, plasma_energy, drude_gamma, oscillators
    )
    constants = optical_constants(energy, epsilon)
    slab = normal_incidence_slab_rt(
        energy,
        np.sqrt(epsilon),
        film_thickness_nm * 1e-9,
        incident_index,
        substrate_index,
    )
    conductivity = optical_conductivity_s_m(energy, epsilon, eps_inf)
    kk_eps1 = kramers_kronig_eps1_finite_window(energy, np.imag(epsilon), eps_inf)
except (ValueError, FloatingPointError) as exc:
    st.error(f"Optical model could not be evaluated: {exc}")
    st.stop()

energy_col, wavelength_col = st.columns(2)
with energy_col:
    st.metric("Spectral window", f"{energy_min:.2f}–{energy_max:.2f} eV")
with wavelength_col:
    st.metric("Equivalent wavelength window", f"{constants['wavelength_m'][-1] * 1e9:.1f}–{constants['wavelength_m'][0] * 1e9:.1f} nm")

tabs = st.tabs(
    [
        "Dielectric function",
        "Optical constants",
        "Thin-film R / T / A",
        "Optical conductivity",
        "Kramers–Kronig check",
    ]
)

with tabs[0]:
    st.subheader("Complex dielectric function")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=energy, y=np.real(epsilon), name="ε₁ (real)", line=dict(width=2.5)))
    fig.add_trace(go.Scatter(x=energy, y=np.imag(epsilon), name="ε₂ (imaginary)", line=dict(width=2.5)))
    fig.add_hline(y=0.0, line_dash="dot")
    fig.update_layout(
        xaxis_title="Photon energy (eV)",
        yaxis_title="Relative dielectric function",
        template="plotly_dark",
        height=500,
    )
    st.plotly_chart(fig, use_container_width=True)
    st.latex(
        r"\epsilon(E)=\epsilon_\infty-\frac{E_p^2}{E^2+i\Gamma_D E}"
        r"+\sum_j\frac{f_j E_p^2}{E_{0j}^2-E^2-i\Gamma_j E}"
    )
    st.caption(
        "The model uses the exp(−iωt) time convention. With non-negative "
        "oscillator strengths and positive damping, ε₂ should be non-negative "
        "for positive photon energies."
    )

with tabs[1]:
    st.subheader("Refractive index and absorption")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=energy, y=constants["n"], name="n", line=dict(width=2.5)))
    fig.add_trace(go.Scatter(x=energy, y=constants["k"], name="k (extinction)", line=dict(width=2.5)))
    fig.update_layout(
        xaxis_title="Photon energy (eV)",
        yaxis_title="Optical constant",
        template="plotly_dark",
        height=440,
    )
    st.plotly_chart(fig, use_container_width=True)

    alpha_cm = constants["absorption_per_m"] / 100.0
    alpha_fig = go.Figure(go.Scatter(x=energy, y=alpha_cm, mode="lines", name="α"))
    alpha_fig.update_layout(
        title="Absorption coefficient",
        xaxis_title="Photon energy (eV)",
        yaxis_title="α (cm⁻¹)",
        template="plotly_dark",
        height=400,
    )
    st.plotly_chart(alpha_fig, use_container_width=True)
    st.latex(r"\tilde{n}=n+ik=\sqrt{\epsilon_1+i\epsilon_2},\qquad "
             r"\alpha=\frac{4\pi k}{\lambda_0}")

with tabs[2]:
    st.subheader("Coherent thin-film optical stack")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=energy, y=slab["reflectance"] * 100.0, name="Reflectance R"))
    fig.add_trace(go.Scatter(x=energy, y=slab["transmittance"] * 100.0, name="Transmittance T"))
    fig.add_trace(go.Scatter(x=energy, y=slab["absorptance"] * 100.0, name="Absorptance A"))
    fig.update_layout(
        xaxis_title="Photon energy (eV)",
        yaxis_title="Power fraction (%)",
        template="plotly_dark",
        height=470,
        yaxis=dict(range=[0, 100]),
    )
    st.plotly_chart(fig, use_container_width=True)
    balance_error = float(np.max(np.abs(
        slab["reflectance"] + slab["transmittance"] + slab["absorptance"] - 1.0
    )))
    st.metric("Maximum energy-balance residual", f"{balance_error:.2e}")
    st.caption(
        "Normal-incidence coherent single-film calculation between lossless "
        "incident and substrate media. The substrate is semi-infinite: its "
        "back interface, roughness, anisotropy and incoherent substrate effects "
        "are not included."
    )

with tabs[3]:
    st.subheader("Complex optical conductivity")
    sigma_real = np.real(conductivity)
    sigma_imag = np.imag(conductivity)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=energy, y=sigma_real / 1e6, name="Re σ", line=dict(width=2.5)))
    fig.add_trace(go.Scatter(x=energy, y=sigma_imag / 1e6, name="Im σ", line=dict(width=2.5)))
    fig.update_layout(
        xaxis_title="Photon energy (eV)",
        yaxis_title="Optical conductivity (10⁶ S m⁻¹)",
        template="plotly_dark",
        height=470,
    )
    st.plotly_chart(fig, use_container_width=True)
    st.latex(r"\sigma(\omega)=-i\epsilon_0\omega"
             r"\left[\epsilon(\omega)-\epsilon_\infty\right]")

with tabs[4]:
    st.subheader("Kramers–Kronig consistency diagnostic")
    valid = np.isfinite(kk_eps1)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=energy, y=np.real(epsilon), name="Drude–Lorentz ε₁", line=dict(width=2.5)
    ))
    fig.add_trace(go.Scatter(
        x=energy[valid], y=kk_eps1[valid], name="Finite-window KK estimate",
        line=dict(width=2, dash="dash")
    ))
    fig.update_layout(
        xaxis_title="Photon energy (eV)",
        yaxis_title="Real dielectric function ε₁",
        template="plotly_dark",
        height=470,
    )
    st.plotly_chart(fig, use_container_width=True)
    st.warning(
        "The KK curve is a finite-window numerical estimate. Missing low- and "
        "high-energy spectral tails can cause disagreement, especially near the "
        "window boundaries. This is a diagnostic, not a proof of causality."
    )

data = pd.DataFrame(
    {
        "photon_energy_eV": energy,
        "vacuum_wavelength_nm": constants["wavelength_m"] * 1e9,
        "epsilon1": np.real(epsilon),
        "epsilon2": np.imag(epsilon),
        "n": constants["n"],
        "k": constants["k"],
        "absorption_cm-1": constants["absorption_per_m"] / 100.0,
        "reflectance": slab["reflectance"],
        "transmittance": slab["transmittance"],
        "absorptance": slab["absorptance"],
        "sigma1_S_m": np.real(conductivity),
        "sigma2_S_m": np.imag(conductivity),
        "epsilon1_KK_finite_window": kk_eps1,
    }
)
st.download_button(
    "Download complete optical spectrum (CSV)",
    data=data.to_csv(index=False).encode("utf-8"),
    file_name="physicsx_optical_spectroscopy.csv",
    mime="text/csv",
)

with st.expander("Scientific scope and assumptions"):
    st.markdown(
        """
        - Linear, local, isotropic dielectric response.
        - Drude free-carrier term plus user-defined Lorentz bound-electron oscillators.
        - Positive photon-energy grid; all model parameters are expressed in eV.
        - The refractive-index branch is selected so n and the extinction coefficient k are non-negative for passive media.
        - Film optics use a coherent plane-parallel layer at normal incidence, between lossless media.
        - No anisotropy, excitons beyond the Lorentz oscillator representation, spatial dispersion, temperature dependence, surface roughness, or realistic material-specific oscillator fit is supplied by default.
        - Results are parameterized model outputs, not experimental data or a claim about a particular material.
        """
    )

st.caption("PHYSICS-X v0.4.0 · Optical spectroscopy and dielectric-response laboratory")
