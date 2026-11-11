# PHYSICS-X

### An Open-Source Computational Physics & Intelligent Research Laboratory

**PHYSICS-X** is a modular, Streamlit-based scientific computing workspace for exploring physical models, running reproducible calculations, and visualizing results.

> Built by **Aman Kumar Patel** with ❤️

---

## Project status

**Early development — v0.1.0**

This first milestone establishes the application architecture and three interactive teaching/research tools:

- **Mechanics Lab:** axial stress, engineering strain, elongation, factor of safety, and an explicitly idealized stress–strain visualization.
- **Quantum Lab:** energy spectrum and stationary wavefunctions of a one-dimensional infinite square well.
- **Crystal Lab:** simple cubic, body-centred cubic, and face-centred cubic unit-cell coordinates, packing fraction, and theoretical density.

These modules are educational and exploratory tools. They are not substitutes for validated engineering design software, experimental measurements, or peer-reviewed numerical workflows.

## Run locally

Requires Python 3.10 or newer.

```bash
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
```

Streamlit will print a local URL (usually http://localhost:8501).

## Repository structure

```text
PHYSICS-X/
├── app.py
├── pages/
│   ├── 1_Mechanics_Lab.py
│   ├── 2_Quantum_Lab.py
│   └── 3_Crystal_Lab.py
├── physicsx/
│   ├── __init__.py
│   ├── mechanics.py
│   ├── quantum.py
│   └── crystal.py
├── .streamlit/
│   └── config.toml
├── .gitignore
└── requirements.txt
```

## Scientific conventions

- SI units are used internally wherever practical.
- Inputs and outputs are labelled with units.
- Simplifying assumptions are stated alongside each model.
- The stress–strain curve in Mechanics Lab is a configurable illustrative bilinear engineering model, not a universal material law.
- The Quantum Lab uses the analytic non-relativistic infinite-well solution.
- The Crystal Lab assumes ideal bulk unit cells and uses the conventional-cell definitions for SC, BCC, and FCC.

## Roadmap

- [x] Streamlit application shell and navigation
- [x] Mechanics: axial loading and idealized stress–strain model
- [x] Quantum: 1D infinite square well
- [x] Crystallography: SC, BCC, FCC unit cells
- [ ] Automated tests and continuous integration
- [ ] Materials data import and provenance tracking
- [ ] Spectroscopy and lattice-dynamics tools
- [ ] Nanophotonics simulation interfaces
- [ ] Reproducible report export
- [ ] Carefully validated machine-learning workflows

## Contributing

Issues, reproducible bug reports, scientific corrections, and pull requests are welcome. Please include units, assumptions, and references for scientific changes.

## License

This project is intended to be released under the Apache License 2.0. A license file will be added with the first tagged release.

---

**PHYSICS-X** — explore the model, inspect the assumptions, reproduce the result.
