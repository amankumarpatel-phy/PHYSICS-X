from __future__ import annotations

import streamlit as st

st.set_page_config(
    page_title="PHYSICS-X | Computational Physics Laboratory",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("⚛️ PHYSICS-X")
st.subheader("Computational Physics & Intelligent Research Laboratory")
st.caption("An open-source scientific workspace • Built by Aman Kumar Patel with ❤️")

st.markdown(
    """
    PHYSICS-X brings interactive physical models and scientific visualizations
    into one modular workspace. Select a laboratory from the sidebar to begin.
    """
)

st.divider()
left, middle, right = st.columns(3)
with left:
    st.metric("Interactive laboratories", "4")
    st.caption("Mechanics • Quantum • Crystallography • Phonons")
with middle:
    st.metric("Release", "0.2.0")
    st.caption("Lattice-dynamics milestone")
with right:
    st.metric("Core language", "Python")
    st.caption("Streamlit interface")

st.header("Available laboratories")
cards = [
    ("Mechanics Lab", "Axial stress, strain, elongation, safety factor and an idealized stress–strain curve."),
    ("Quantum Lab", "Energy levels and stationary states for the one-dimensional infinite square well."),
    ("Crystal Lab", "SC, BCC and FCC unit cells, packing fraction, atomic radius and theoretical density."),
    ("Phonon Lab", "Acoustic and optical dispersion branches plus a sampled 1D phonon density of states."),
]
for title, description in cards:
    with st.container(border=True):
        st.markdown(f"### {title}")
        st.write(description)

st.info(
    "Research integrity: every model has a stated scope and assumptions. "
    "Check the model description before interpreting results."
)

with st.expander("Project roadmap"):
    st.markdown(
        """
        - [x] Streamlit application shell
        - [x] Mechanics, quantum-well and cubic-crystal models
        - [x] 1D diatomic-chain phonon dispersion and density of states
        - [x] Unit tests and GitHub Actions continuous integration
        - [ ] 3D lattice dynamics and force-constant input
        - [ ] Spectroscopy and data provenance
        - [ ] Nanophotonics simulation workflows
        - [ ] Reproducible report export
        - [ ] Validated AI-assisted optimization
        """
    )

st.markdown("---")
st.caption("PHYSICS-X v0.2.0 · Open-source computational physics project")
