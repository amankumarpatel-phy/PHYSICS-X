from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from physicsx.crystal import STRUCTURES, atomic_radius_m, theoretical_density_kg_m3, unit_cell_points

st.set_page_config(page_title="Crystal Lab | PHYSICS-X", page_icon="💎", layout="wide")
st.title("💎 Crystal Lab")
st.write("Inspect ideal cubic unit cells and calculate basic geometric properties.")

structure = st.selectbox("Crystal structure", list(STRUCTURES.keys()))
col1, col2 = st.columns(2)
with col1:
    lattice_angstrom = st.number_input("Lattice constant (Å)", min_value=0.01, value=3.61, step=0.01)
with col2:
    molar_mass_g_mol = st.number_input("Molar mass (g mol⁻¹)", min_value=0.001, value=63.546, step=0.1)

info = STRUCTURES[structure]
a_m = lattice_angstrom * 1e-10
points = unit_cell_points(structure, a_m)
radius = atomic_radius_m(structure, a_m)
density = theoretical_density_kg_m3(structure, a_m, molar_mass_g_mol * 1e-3)

m1, m2, m3, m4 = st.columns(4)
m1.metric("Atoms / conventional cell", str(info["atoms_per_cell"]))
m2.metric("Coordination number", str(info["coordination_number"]))
m3.metric("Atomic radius", f"{radius * 1e10:.4f} Å")
m4.metric("Packing fraction", f"{info['packing_fraction'] * 100:.2f}%")

st.metric("Theoretical density", f"{density / 1000:.4f} g cm⁻³")

# Draw the unit-cell edges.
corners = np.array(
    [[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)],
    dtype=float,
) * a_m
edges = []
for i, p in enumerate(corners):
    for q in corners[i + 1:]:
        if np.sum(np.abs(p - q) > 0) == 1:
            edges.append((p, q))

fig = go.Figure()
for p, q in edges:
    fig.add_trace(go.Scatter3d(
        x=[p[0] / a_m, q[0] / a_m],
        y=[p[1] / a_m, q[1] / a_m],
        z=[p[2] / a_m, q[2] / a_m],
        mode="lines",
        line=dict(width=5),
        showlegend=False,
        hoverinfo="skip",
    ))
fig.add_trace(go.Scatter3d(
    x=points[:, 0] / a_m,
    y=points[:, 1] / a_m,
    z=points[:, 2] / a_m,
    mode="markers",
    marker=dict(size=9),
    name="Lattice points",
    text=[f"({x:.1f}, {y:.1f}, {z:.1f})" for x, y, z in points / a_m],
    hovertemplate="Fractional coordinate: %{text}<extra></extra>",
))
fig.update_layout(
    title=f"{structure} conventional unit cell",
    scene=dict(
        xaxis_title="x / a",
        yaxis_title="y / a",
        zaxis_title="z / a",
        aspectmode="cube",
    ),
    template="plotly_dark",
    height=600,
)
st.plotly_chart(fig, use_container_width=True)

with st.expander("Crystallographic relations"):
    st.latex(r"\rho = \frac{Z M}{N_A a^3}")
    st.latex(r"\text{Packing fraction} = \frac{Z(4\pi r^3/3)}{a^3}")
    st.markdown(
        """
        Density uses the conventional-cell atom count Z and assumes a pure,
        fully occupied, ideal bulk crystal. Packing fraction and radius use the
        hard-sphere contact geometry of the selected cubic lattice. Real crystals
        may contain defects, vacancies, thermal expansion, multiple species or
        non-ideal occupancies; those effects are not represented here.
        """
    )
