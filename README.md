# PHYSICS-X

### An Open-Source Computational Physics & Intelligent Research Laboratory

**PHYSICS-X** is a modular, Streamlit-based scientific computing workspace for exploring physical models, running reproducible calculations, and visualizing results.

> Built by **Aman Kumar Patel** with ❤️

---

## Project status

**Early development — v0.4.0**

The current release contains six interactive laboratories:

- **Mechanics Lab:** axial stress, engineering strain, elongation, factor of safety, and an explicitly idealized stress–strain visualization.
- **Quantum Lab:** energy spectrum and stationary wavefunctions of a one-dimensional infinite square well.
- **Crystal Lab:** simple cubic, body-centred cubic, and face-centred cubic unit-cell coordinates, packing fraction, and theoretical density.
- **Phonon Lab:** acoustic and optical dispersion branches and a sampled phonon density of states for a one-dimensional diatomic chain.
- **Reciprocal Space Lab:** SC, BCC, FCC and HCP geometry; Miller-plane sections; reciprocal-lattice points; a Wigner–Seitz construction of the first Brillouin zone; and a simplified powder XRD pattern.
- **Optical Spectroscopy Lab:** Drude–Lorentz dielectric response, complex refractive index, absorption coefficient, coherent thin-film R/T/A, optical conductivity and a finite-window Kramers–Kronig diagnostic.

The optical module is a parameterized model, not an electronic-structure calculation or experimental material database. The XRD engine assumes equal atomic form factors and Gaussian peak broadening; it is intended for teaching and algorithm development, not quantitative phase identification.

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
│   ├── 5_Reciprocal_Space_Lab.py
│   └── 6_Optical_Spectroscopy_Lab.py
├── physicsx/
│   ├── __init__.py
│   ├── mechanics.py
│   ├── quantum.py
│   ├── crystal.py
│   ├── crystal_advanced.py
│   ├── phonon.py
│   └── optics.py
├── tests/
│   ├── test_models.py
│   ├── test_phonon.py
│   ├── test_crystal_advanced.py
│   └── test_optics.py
├── .github/workflows/
│   └── tests.yml
├── .streamlit/
│   └── config.toml
├── LICENSE
├── .gitignore
└── requirements.txt
~~~

## Scientific conventions

- SI units are used internally for optical conductivity, wavelength, absorption and layer thickness; UI inputs are clearly labelled.
- The optical dielectric model is expressed in photon-energy units E = ℏω and uses the exp(−iωt) convention:
  ε(E) = ε∞ − Ep²/(E² + iΓD E) + Σ fⱼEp²/(E₀ⱼ² − E² − iΓⱼE).
- Passive optical constants use the refractive-index branch n ≥ 0 and k ≥ 0. The absorption coefficient is α = 4πk/λ₀.
- Thin-film R/T/A is a coherent, normal-incidence, single-layer transfer-matrix calculation between lossless incident and substrate media. The substrate is semi-infinite.
- Optical conductivity is calculated as σ(ω) = −iε₀ω[ε(ω) − ε∞].
- The Kramers–Kronig display is a finite-window numerical principal-value estimate. Spectral tails outside the selected energy interval are omitted; it is a diagnostic rather than an exact infinite-range reconstruction.
- The Phonon Lab uses a harmonic one-dimensional diatomic chain with two atoms per primitive cell and equal nearest-neighbour springs.
- The Reciprocal Space Lab uses row-wise direct and reciprocal basis vectors satisfying aᵢ·bⱼ = 2πδᵢⱼ. Its XRD pattern uses equal atomic form factors and Gaussian broadening.

## Roadmap

- [x] Streamlit application shell and navigation
- [x] Mechanics: axial loading and idealized stress–strain model
- [x] Quantum: 1D infinite square well
- [x] Crystallography: SC, BCC, FCC and HCP geometry
- [x] 1D diatomic-chain phonon dispersion and DOS
- [x] Miller-plane sections and reciprocal-lattice visualization
- [x] First Brillouin-zone Wigner–Seitz construction
- [x] Simplified powder XRD pattern generation
- [x] Drude–Lorentz optical spectra and complex optical constants
- [x] Coherent thin-film R/T/A and complex optical conductivity
- [x] Finite-window Kramers–Kronig diagnostic
- [x] Unit tests and GitHub Actions continuous integration
- [ ] Full 3D lattice dynamics with realistic force constants
- [ ] Atomic form-factor library and quantitative diffraction refinements
- [ ] Measured-spectrum fitting and material parameter provenance
- [ ] Nanophotonics simulation interfaces
- [ ] Reproducible report export
- [ ] Carefully validated machine-learning workflows

## Contributing

Issues, reproducible bug reports, scientific corrections, and pull requests are welcome. Please include units, assumptions, and references for scientific changes.

## License

Released under the Apache License 2.0. See [LICENSE](LICENSE) for details.

---

**PHYSICS-X** — explore the model, inspect the assumptions, reproduce the result.
