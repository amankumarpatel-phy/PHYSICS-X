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
    st.metric("Interactive laboratories", "3")
    st.caption("Mechanics • Quantum • Crystallography")
with middle:
    st.metric("Release", "0.1.0")
    st.caption("Early development")
with right:
    st.metric("Core language", "Python")
    st.caption("Streamlit interface")

st.header("Available laboratories")
cards = [
    ("Mechanics Lab", "Axial stress, strain, elongation, safety factor and an idealized stress–strain curve."),
    ("Quantum Lab", "Energy levels and stationary states for the one-dimensional infinite square well."),
    ("Crystal Lab", "SC, BCC and FCC unit cells, packing fraction, atomic radius and theoretical density."),
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
        - [x] Elementary mechanics, quantum-well and cubic-crystal models
        - [ ] Automated scientific tests
        - [ ] Spectroscopy and lattice dynamics
        - [ ] Nanophotonics simulation workflows
        - [ ] Reproducible export and data provenance
        - [ ] Validated AI-assisted optimization
        """
    )

st.markdown("---")
st.caption("PHYSICS-X v0.1.0 · Open-source computational physics project")
