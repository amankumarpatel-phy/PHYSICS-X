from __future__ import annotations

from itertools import product

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from physicsx.crystal_advanced import (
    first_brillouin_zone,
    fractional_to_cartesian,
    lattice_geometry,
    miller_d_spacing_m,
    miller_plane_section,
    powder_xrd_pattern,
    reciprocal_lattice_points,
)

st.set_page_config(
    page_title="Reciprocal Space Lab | PHYSICS-X",
    page_icon="🧭",
    layout="wide",
)
st.title("🧭 3D Crystallography & Reciprocal Space")
st.write(
    "Explore unit-cell geometry, Miller planes, reciprocal-lattice points, "
    "the first Brillouin zone and a simplified powder X-ray diffraction pattern."
)

with st.sidebar:
    st.header("Crystal parameters")
    structure = st.selectbox(
        "Bravais lattice / structure",
        [
            "Simple Cubic (SC)",
            "Body-Centred Cubic (BCC)",
            "Face-Centred Cubic (FCC)",
            "Hexagonal Close-Packed (HCP)",
        ],
    )
    a_angstrom = st.number_input(
        "Lattice parameter a (Å)", min_value=0.1, max_value=30.0, value=3.61, step=0.01
    )
    if structure == "Hexagonal Close-Packed (HCP)":
        c_over_a = st.number_input(
            "Axial ratio c/a", min_value=0.5, max_value=3.0, value=1.633, step=0.001
        )
        c_angstrom = a_angstrom * c_over_a
    else:
        c_angstrom = None
        c_over_a = None

try:
    geometry = lattice_geometry(
        structure,
        a_angstrom * 1e-10,
        None if c_angstrom is None else c_angstrom * 1e-10,
    )
except ValueError as exc:
    st.error(str(exc))
    st.stop()

direct = geometry["direct"]
reciprocal = geometry["reciprocal"]
basis = geometry["basis"]
cell_volume = float(geometry["cell_volume_m3"])
primitive_volume = float(geometry["primitive_volume_m3"])

st.caption(
    f"Cell volume: {cell_volume * 1e30:.4f} Å³ · "
    f"Atoms in displayed cell: {geometry['atoms_per_cell']} · "
    f"Primitive volume: {primitive_volume * 1e30:.4f} Å³"
)

tab_cell, tab_recip, tab_bz, tab_xrd = st.tabs(
    ["Unit cell & Miller planes", "Reciprocal lattice", "First Brillouin zone", "Powder XRD"]
)

with tab_cell:
    st.subheader("Interactive unit cell")
    hcol, kcol, lcol = st.columns(3)
    with hcol:
        h = st.number_input("Miller index h", min_value=-8, max_value=8, value=1, step=1)
    with kcol:
        k = st.number_input("Miller index k", min_value=-8, max_value=8, value=1, step=1)
    with lcol:
        ell = st.number_input("Miller index l", min_value=-8, max_value=8, value=1, step=1)

    hkl = (int(h), int(k), int(ell))
    if hkl == (0, 0, 0):
        st.error("(000) is not a valid Miller plane.")
    else:
        plane_values = np.array(list(product((0.0, 1.0), repeat=3))) @ np.asarray(hkl, dtype=float)
        plane_min, plane_max = float(np.min(plane_values)), float(np.max(plane_values))
        default_offset = min(max(0.5 * (plane_min + plane_max), plane_min + 1e-6), plane_max - 1e-6)
        offset = st.slider(
            "Plane offset: h·u + k·v + l·w",
            min_value=float(plane_min),
            max_value=float(plane_max),
            value=float(default_offset),
        )
        try:
            section = miller_plane_section(hkl, offset, direct)
            section_frac = section @ np.linalg.inv(direct)
            corners_frac = np.array(list(product((0.0, 1.0), repeat=3)))
            corners = fractional_to_cartesian(corners_frac, direct)
            atoms = fractional_to_cartesian(basis, direct)

            fig = go.Figure()
            # Twelve edges of the fractional parallelepiped.
            for i, first in enumerate(corners_frac):
                for second in corners_frac[i + 1:]:
                    if np.count_nonzero(np.abs(first - second) > 0.0) == 1:
                        edge = fractional_to_cartesian(np.vstack([first, second]), direct)
                        fig.add_trace(go.Scatter3d(
                            x=edge[:, 0] * 1e10, y=edge[:, 1] * 1e10, z=edge[:, 2] * 1e10,
                            mode="lines", line=dict(width=5), showlegend=False, hoverinfo="skip",
                        ))
            count = len(section)
            fig.add_trace(go.Mesh3d(
                x=section[:, 0] * 1e10,
                y=section[:, 1] * 1e10,
                z=section[:, 2] * 1e10,
                i=[0] * (count - 2),
                j=list(range(1, count - 1)),
                k=list(range(2, count)),
                opacity=0.55,
                color="#35C8B2",
                name=f"({h} {k} {ell}) plane",
                hovertemplate="Plane section<extra></extra>",
            ))
            fig.add_trace(go.Scatter3d(
                x=atoms[:, 0] * 1e10, y=atoms[:, 1] * 1e10, z=atoms[:, 2] * 1e10,
                mode="markers", marker=dict(size=9, color="#F4B942"),
                name="Basis atoms",
            ))
            fig.update_layout(
                title=f"{structure}: fractional plane ({h} {k} {ell})",
                scene=dict(xaxis_title="x (Å)", yaxis_title="y (Å)", zaxis_title="z (Å)", aspectmode="data"),
                template="plotly_dark", height=620,
            )
            st.plotly_chart(fig, use_container_width=True)
            d_hkl = miller_d_spacing_m(hkl, reciprocal)
            st.metric("Interplanar spacing d(hkl)", f"{d_hkl * 1e10:.5f} Å")
            st.caption(
                "The plane is a geometric section of the displayed cell. For centred lattices, "
                "not every geometrically indexed plane corresponds to an allowed diffraction reflection."
            )
        except ValueError as exc:
            st.warning(str(exc))

with tab_recip:
    st.subheader("Reciprocal-lattice vectors")
    st.latex(r"\mathbf{b}_i\cdot\mathbf{a}_j=2\pi\delta_{ij}")
    labels = ["b₁", "b₂", "b₃"]
    vector_rows = []
    for index, vector in enumerate(np.asarray(reciprocal)):
        vector_rows.append({
            "Vector": labels[index],
            "x (Å⁻¹)": vector[0] * 1e-10,
            "y (Å⁻¹)": vector[1] * 1e-10,
            "z (Å⁻¹)": vector[2] * 1e-10,
            "Magnitude (Å⁻¹)": np.linalg.norm(vector) * 1e-10,
        })
    st.dataframe(pd.DataFrame(vector_rows).set_index("Vector"), use_container_width=True)

    indices, vectors = reciprocal_lattice_points(
        np.asarray(reciprocal), structure, np.asarray(basis), index_limit=2
    )
    fig_recip = go.Figure()
    fig_recip.add_trace(go.Scatter3d(
        x=vectors[:, 0] * 1e-10,
        y=vectors[:, 1] * 1e-10,
        z=vectors[:, 2] * 1e-10,
        mode="markers",
        marker=dict(size=4),
        text=[f"({i} {j} {k})" for i, j, k in indices],
        hovertemplate="Reciprocal point %{text}<extra></extra>",
        name="Allowed reciprocal points",
    ))
    fig_recip.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0], mode="markers",
        marker=dict(size=7, symbol="diamond", color="#F4B942"), name="Γ",
    ))
    fig_recip.update_layout(
        title="Reciprocal lattice (2π convention)",
        scene=dict(xaxis_title="kₓ (Å⁻¹)", yaxis_title="kᵧ (Å⁻¹)", zaxis_title="k_z (Å⁻¹)", aspectmode="data"),
        template="plotly_dark", height=600,
    )
    st.plotly_chart(fig_recip, use_container_width=True)
    st.caption(
        "Cubic BCC/FCC reciprocal points are filtered using the centring selection rules. "
        "HCP is treated as a primitive hexagonal Bravais lattice with a two-atom basis."
    )

with tab_bz:
    st.subheader("First Brillouin zone")
    st.write("The first Brillouin zone is the Wigner–Seitz cell of the reciprocal lattice.")
    try:
        bz = first_brillouin_zone(
            np.asarray(reciprocal), structure, np.asarray(basis), index_limit=3
        )
        vertices = np.asarray(bz["vertices"])
        triangles = np.asarray(bz["simplices"])
        fig_bz = go.Figure()
        fig_bz.add_trace(go.Mesh3d(
            x=vertices[:, 0] * 1e-10,
            y=vertices[:, 1] * 1e-10,
            z=vertices[:, 2] * 1e-10,
            i=triangles[:, 0],
            j=triangles[:, 1],
            k=triangles[:, 2],
            opacity=0.38,
            color="#3D9BE9",
            name="First Brillouin zone",
            hovertemplate="BZ boundary vertex<extra></extra>",
        ))
        fig_bz.add_trace(go.Scatter3d(
            x=vertices[:, 0] * 1e-10,
            y=vertices[:, 1] * 1e-10,
            z=vertices[:, 2] * 1e-10,
            mode="markers",
            marker=dict(size=3, color="#F4B942"),
            name="Vertices",
        ))
        fig_bz.add_trace(go.Scatter3d(
            x=[0], y=[0], z=[0],
            mode="markers",
            marker=dict(size=7, color="#F4B942", symbol="diamond"),
            name="Γ",
        ))
        fig_bz.update_layout(
            title=f"First Brillouin zone · {structure}",
            scene=dict(xaxis_title="kₓ (Å⁻¹)", yaxis_title="kᵧ (Å⁻¹)", zaxis_title="k_z (Å⁻¹)", aspectmode="data"),
            template="plotly_dark", height=620,
        )
        st.plotly_chart(fig_bz, use_container_width=True)

        expected_bz_volume = (2.0 * np.pi) ** 3 / primitive_volume
        relative_error = abs(float(bz["volume_m3"]) - expected_bz_volume) / expected_bz_volume
        c1, c2 = st.columns(2)
        c1.metric("Numerical BZ volume", f"{float(bz['volume_m3']) * 1e30:.5f} Å⁻³")
        c2.metric("Volume identity relative error", f"{relative_error:.2e}")
        st.caption(
            "The volume identity is V_BZ = (2π)³ / V_primitive. The zone is constructed "
            "from reciprocal vectors within the selected integer-index shell; the identity "
            "provides a useful numerical completeness check."
        )
    except Exception as exc:
        st.error(f"Brillouin-zone construction failed: {exc}")

with tab_xrd:
    st.subheader("Powder X-ray diffraction")
    xcol, ycol = st.columns(2)
    with xcol:
        wavelength_angstrom = st.selectbox(
            "X-ray wavelength (Å)",
            [("Cu Kα₁ (1.5406 Å)", 1.5406), ("Mo Kα₁ (0.7093 Å)", 0.7093)],
            format_func=lambda item: item[0],
        )[1]
    with ycol:
        fwhm = st.slider("Gaussian peak FWHM (° 2θ)", 0.05, 0.8, 0.16, 0.01)
    try:
        angle, intensity, peaks = powder_xrd_pattern(
            np.asarray(reciprocal),
            np.asarray(basis),
            structure,
            wavelength_m=wavelength_angstrom * 1e-10,
            max_index=7,
            fwhm_deg=fwhm,
        )
        fig_xrd = go.Figure()
        fig_xrd.add_trace(go.Scatter(
            x=angle, y=intensity, mode="lines", name="Simulated pattern", line=dict(width=2)
        ))
        fig_xrd.update_layout(
            title="Simulated powder diffraction pattern",
            xaxis_title="Diffraction angle 2θ (degrees)",
            yaxis_title="Relative intensity (a.u.)",
            template="plotly_dark", height=480,
        )
        st.plotly_chart(fig_xrd, use_container_width=True)

        if peaks:
            peak_table = pd.DataFrame(peaks)
            peak_table["intensity"] = peak_table["intensity"].astype(float)
            peak_table = peak_table.sort_values("intensity", ascending=False).head(20)
            st.subheader("Strongest simulated peaks")
            st.dataframe(peak_table, use_container_width=True, hide_index=True)
            st.download_button(
                "Download peak list (CSV)",
                data=pd.DataFrame(peaks).to_csv(index=False).encode("utf-8"),
                file_name="physicsx_powder_xrd_peaks.csv",
                mime="text/csv",
            )
        else:
            st.warning("No reflections fell within the selected 2θ range for these parameters.")
    except ValueError as exc:
        st.error(str(exc))

    with st.expander("Diffraction model and limitations"):
        st.latex(r"2d_{hkl}\sin\theta=n\lambda")
        st.latex(r"F_{hkl}=\sum_j f_j\exp\!\left[2\pi i(hx_j+ky_j+lz_j)\right]")
        st.markdown(
            """
            This is a teaching-level, kinematic powder pattern. It uses first-order
            Bragg diffraction, equal atomic form factors (f = 1), structure-factor
            intensities and Gaussian broadening. It omits Lorentz-polarization factors,
            multiplicity corrections beyond the explicitly summed reciprocal points,
            Debye-Waller factors, anomalous scattering, absorption, preferred orientation,
            instrumental profile convolution and realistic atomic scattering factors.
            It should not be used for quantitative phase identification.
            """
        )

st.caption("PHYSICS-X v0.3.0 · Crystallography and reciprocal-space tools")
