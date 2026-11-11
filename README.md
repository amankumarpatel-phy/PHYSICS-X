# PHYSICS-X

### An Open-Source Computational Physics & Intelligent Research Laboratory

**PHYSICS-X** is a modular, Streamlit-based scientific computing workspace for exploring physical models, running reproducible calculations, and visualizing results.

> Built by **Aman Kumar Patel** with ❤️

---

## Project status

**Early development — v0.3.0**

The current release contains five interactive laboratories:

- **Mechanics Lab:** axial stress, engineering strain, elongation, factor of safety, and an explicitly idealized stress–strain visualization.
- **Quantum Lab:** energy spectrum and stationary wavefunctions of a one-dimensional infinite square well.
- **Crystal Lab:** simple cubic, body-centred cubic, and face-centred cubic unit-cell coordinates, packing fraction, and theoretical density.
- **Phonon Lab:** acoustic and optical dispersion branches and a sampled phonon density of states for a one-dimensional diatomic chain.
- **Reciprocal Space Lab:** SC, BCC, FCC and HCP geometry; Miller-plane sections; reciprocal-lattice points; a Wigner–Seitz construction of the first Brillouin zone; and a simplified powder XRD pattern.

The crystallography extension uses a 2π reciprocal-lattice convention. The XRD engine assumes equal atomic form factors and Gaussian peak broadening; it is intended for teaching and algorithm development, not quantitative phase identification.

These modules are educational and exploratory tools. They are not substitutes for validated engineering design software, experimental measurements, or peer-reviewed numerical workflows.

## Run locally

Requires Python 3.10 or newer.

~~~bash
git clone https://github.com/amankumarpatel-phy/PHYSICS-X.git
cd PHYSICS-X

python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

python -m pip install --upgrade pip
pip install -r requirements.txt

streamlit run app.py
~~~

Streamlit will print a local URL (usually http://localhost:8501).

## Repository structure

~~~text
PHYSICS-X/
├── app.py
├── pages/
│   ├── 1_Mechanics_Lab.py
│   ├── 2_Quantum_Lab.py
│   ├── 3_Crystal_Lab.py
│   ├── 4_Phonon_Lab.py
│   └── 5_Reciprocal_Space_Lab.py
├── physicsx/
│   ├── __init__.py
│   ├── mechanics.py
│   ├── quantum.py
│   ├── crystal.py
│   ├── crystal_advanced.py
│   └── phonon.py
├── tests/
│   ├── test_models.py
│   ├── test_phonon.py
│   └── test_crystal_advanced.py
├── .github/workflows/
│   └── tests.yml
├── .streamlit/
│   └── config.toml
├── LICENSE
├── .gitignore
└── requirements.txt
~~~

## Scientific conventions

- SI units are used internally wherever practical.
- Inputs and outputs are labelled with units.
- Simplifying assumptions are stated alongside each model.
- The stress–strain curve in Mechanics Lab is a configurable illustrative bilinear engineering model, not a universal material law.
- The Quantum Lab uses the analytic non-relativistic infinite-well solution.
- The Crystal Lab assumes ideal bulk unit cells and uses conventional-cell definitions for SC, BCC, and FCC.
- The Phonon Lab uses a harmonic, one-dimensional diatomic chain with two atoms per primitive cell and equal nearest-neighbour springs. Its DOS is a sampled histogram normalized to two modes per primitive cell.
- The Reciprocal Space Lab uses row-wise direct and reciprocal basis vectors satisfying aᵢ·bⱼ = 2πδᵢⱼ. Its XRD pattern uses equal atomic form factors, structure-factor intensities, summed reciprocal-grid multiplicities and Gaussian broadening. Corrections for realistic atomic scattering, Debye–Waller effects, absorption and instrument response are not yet included.

## Roadmap

- [x] Streamlit application shell and navigation
- [x] Mechanics: axial loading and idealized stress–strain model
- [x] Quantum: 1D infinite square well
- [x] Crystallography: SC, BCC, FCC and HCP geometry
- [x] 1D diatomic-chain phonon dispersion and DOS
- [x] Miller-plane sections and reciprocal-lattice visualization
- [x] First Brillouin-zone Wigner–Seitz construction
- [x] Simplified powder XRD pattern generation
- [x] Unit tests and GitHub Actions continuous integration
- [ ] Full 3D lattice dynamics with realistic force constants
- [ ] Atomic form-factor library and quantitative diffraction refinements
- [ ] Spectroscopy and optical response tools
- [ ] Nanophotonics simulation interfaces
- [ ] Reproducible report export
- [ ] Carefully validated machine-learning workflows

## Contributing

Issues, reproducible bug reports, scientific corrections, and pull requests are welcome. Please include units, assumptions, and references for scientific changes.

## License

Released under the Apache License 2.0. See [LICENSE](LICENSE) for details.

---

**PHYSICS-X** — explore the model, inspect the assumptions, reproduce the result.
