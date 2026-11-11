from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from physicsx.mechanics import axial_response, idealized_bilinear_curve

st.set_page_config(page_title="Mechanics Lab | PHYSICS-X", page_icon="🏗️", layout="wide")
st.title("🏗️ Mechanics Lab")
st.write("Explore the response of a straight, prismatic bar under uniform axial loading.")

with st.sidebar:
    st.header("Specimen and loading")
    force_kn = st.number_input("Axial force (kN)", min_value=-10000.0, max_value=10000.0, value=25.0, step=1.0)
    area_mm2 = st.number_input("Cross-sectional area (mm²)", min_value=0.01, value=250.0, step=10.0)
    length_mm = st.number_input("Original length (mm)", min_value=0.1, value=1000.0, step=50.0)
    modulus_gpa = st.number_input("Young's modulus (GPa)", min_value=0.001, value=200.0, step=5.0)
    yield_mpa = st.number_input("Yield strength (MPa)", min_value=0.001, value=250.0, step=10.0)
    ultimate_mpa = st.number_input("Ultimate strength (MPa)", min_value=0.001, value=400.0, step=10.0)
    ultimate_strain = st.number_input("Strain at ultimate point", min_value=0.00001, value=0.15, format="%.5f")

try:
    result = axial_response(
        force_n=force_kn * 1e3,
        area_m2=area_mm2 * 1e-6,
        length_m=length_mm * 1e-3,
        youngs_modulus_pa=modulus_gpa * 1e9,
        yield_strength_pa=yield_mpa * 1e6,
    )
    strain, stress, yield_strain = idealized_bilinear_curve(
        modulus_gpa * 1e9, yield_mpa * 1e6, ultimate_mpa * 1e6, ultimate_strain
    )
except ValueError as exc:
    st.error(str(exc))
    st.stop()

c1, c2, c3, c4 = st.columns(4)
c1.metric("Engineering stress", f"{result['stress_pa'] / 1e6:.3f} MPa")
c2.metric("Engineering strain", f"{result['strain']:.4e}")
c3.metric("Axial deformation", f"{result['elongation_m'] * 1e3:.4f} mm")
c4.metric("Factor of safety", "∞" if np.isinf(result["factor_of_safety"]) else f"{result['factor_of_safety']:.3f}")

if result["yield_utilization"] >= 1.0:
    st.warning("The calculated elastic stress reaches or exceeds the specified yield strength. The linear-elastic deformation result should not be treated as a valid post-yield prediction.")
else:
    st.success("The calculated stress is below the specified yield strength under the stated linear-elastic assumptions.")

fig = go.Figure()
fig.add_trace(go.Scatter(x=strain, y=stress / 1e6, mode="lines", name="Illustrative material curve", line=dict(width=3)))
fig.add_vline(x=yield_strain, line_dash="dash", annotation_text="Specified yield point")
fig.add_hline(y=abs(result["stress_pa"]) / 1e6, line_dash="dot", annotation_text="Applied stress magnitude")
fig.update_layout(
    title="Illustrative engineering stress–strain response",
    xaxis_title="Engineering strain (dimensionless)",
    yaxis_title="Engineering stress (MPa)",
    template="plotly_dark",
    height=520,
)
st.plotly_chart(fig, use_container_width=True)

with st.expander("Model assumptions and equations"):
    st.latex(r"\sigma = \frac{F}{A}")
    st.latex(r"\epsilon = \frac{\sigma}{E}")
    st.latex(r"\Delta L = \epsilon L")
    st.latex(r"n = \frac{\sigma_y}{|\sigma|}")
    st.markdown(
        """
        The bar calculation assumes uniform uniaxial stress, a constant cross-section,
        small deformation and linear elasticity. Compression is represented by a
        negative force and stress. The curve is an illustrative bilinear construction:
        elastic up to the specified yield point, followed by a straight strain-hardening
        segment to the user-specified ultimate point. It does **not** model necking,
        fracture, plastic unloading or a material-specific constitutive response.
        """
    )
