# CP2K Learn — A Beginner's Tutorial for CP2K

[English](README.md) | [简体中文](README_zh.md)

---

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![CP2K Version](https://img.shields.io/badge/CP2K-v2025.2-green.svg)](https://www.cp2k.org/)
[![Chapters](https://img.shields.io/badge/chapters-16-orange.svg)](#tutorial-chapters)

> A hands-on, difficulty-graded CP2K tutorial series with 16 chapters, 30+ tested input files, and an interactive web interface. All examples have been validated against CP2K v2025.2.

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Tutorial Chapters](#tutorial-chapters)
- [Project Structure](#project-structure)
- [Web Interface](#web-interface)
- [Contributing](#contributing)
- [License](#license)

## Introduction

**CP2K Learn** is a structured learning path for [CP2K](https://www.cp2k.org/), a powerful open-source quantum chemistry and solid-state physics software package. This project provides:

- **16 progressive chapters** — from basic single-point energy to advanced GW/BSE band structure calculations
- **30+ validated input files** — every `.inp` file has been tested with `cp2k.ssmp` (CP2K v2025.2)
- **Ready-to-use structure files** — `.xyz` coordinate files for water boxes, silicon crystals, alanine dipeptide, and more
- **Interactive web interface** — a single-page application for browsing tutorials with progress tracking

### What You Will Learn

| Level | Chapters | Topics |
|-------|----------|--------|
| ⭐ Beginner | 1–4 | First CP2K run, DFT basics, basis sets, cutoff convergence |
| ⭐⭐ Intermediate | 5–9 | Periodic systems, geometry optimization, cell optimization, molecular dynamics |
| ⭐⭐⭐ Advanced | 10–14 | Metadynamics, xTB, QM/MM, vibrational analysis, TDDFT |
| ⭐⭐⭐⭐⭐ Expert | 15–16 | Post-Hartree-Fock (PBE0 hybrid), band structure |

## Prerequisites

| Dependency | Version | Required |
|------------|---------|----------|
| CP2K | >= 2025.2 | Yes |
| Python | >= 3.8 | No (for web interface scripts) |

### Installing CP2K

Download the pre-compiled binary from the [CP2K GitHub Releases](https://github.com/cp2k/cp2k/releases) page:

```bash
# Example: download CP2K 2025.2 for Linux x86_64
wget https://github.com/cp2k/cp2k/releases/download/v2025.2/cp2k-2025.2-Linux-gnu-x86_64.ssmp
chmod +x cp2k-2025.2-Linux-gnu-x86_64.ssmp
```

Or install via package manager:

```bash
# conda (recommended for beginners)
conda install -c conda-forge cp2k

# Ubuntu/Debian
sudo apt install cp2k
```

### Data Files

CP2K requires data files for basis sets and pseudopotentials. Set the following in your input files:

```
BASIS_SET_FILE_NAME  BASIS_SET
POTENTIAL_FILE_NAME  GTH_POTENTIALS
```

These files are located in `<CP2K_INSTALL_DIR>/data/`. For MOLOPT basis sets, use `BASIS_MOLOPT`.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/user/cp2k-learn.git
cd cp2k-learn

# Run your first CP2K calculation
cd 01-hello-cp2k
cp2k.ssmp -i helium.inp -o helium.out

# Check the output
grep "ENERGY" helium.out
```

## Tutorial Chapters

### ⭐ Beginner

| # | Chapter | Input Files | Description |
|---|---------|-------------|-------------|
| 01 | [Hello CP2K](01-hello-cp2k/) | `helium.inp` | Your first CP2K calculation — a single helium atom |
| 02 | [Water Static](02-water-static/) | `H2O_static.inp` | Static DFT calculation on a water molecule |
| 03 | [Basis Sets](03-basis-sets/) | `H2O_szv.inp`, `H2O_dzvp.inp`, `H2O_tzv2p.inp` | Comparing SZV, DZVP, and TZVP basis sets |
| 04 | [Cutoff Convergence](04-cutoff-convergence/) | `Si_cutoff.inp` | Convergence testing of plane-wave cutoff energy |

### ⭐⭐ Intermediate

| # | Chapter | Input Files | Description |
|---|---------|-------------|-------------|
| 05 | [Periodic Silicon](05-periodic-silicon/) | `Si_bulk8.inp`, `Si_bulk8_metal.inp` | Periodic boundary conditions, metallic systems |
| 06 | [Geometry Optimization](06-geo-opt/) | `H2O_geo_opt.inp`, `H2O_geo_opt_fixed.inp` | Optimizing atomic positions with CG and BFGS |
| 07 | [Cell Optimization](07-cell-opt/) | `Si_cell_opt.inp` | Optimizing lattice parameters for silicon |
| 08 | [NVE Molecular Dynamics](08-md-nve/) | `H2O_nve.inp` | Microcanonical ensemble MD simulation |
| 09 | [NVT/NPT Thermostats](09-md-thermostats/) | `H2O_nvt.inp`, `H2O_nvt_csvr.inp`, `H2O_npt.inp` | Canonical and isothermal-isobaric MD with CSVR thermostat |

### ⭐⭐⭐ Advanced

| # | Chapter | Input Files | Description |
|---|---------|-------------|-------------|
| 10 | [Metadynamics](10-metadynamics/) | `alanine_dipeptide_meta.inp` | Enhanced sampling with Well-Tempered Metadynamics |
| 11 | [xTB Method](11-xtb/) | `ethanol_xtb.inp`, `water_xtb_md.inp` | GFN1-xTB semi-empirical method for fast calculations |
| 12 | [QM/MM Basics](12-qmmm-basics/) | `water_qmmm.inp` | Hybrid QM/MM simulation of water |
| 13 | [Vibrational Analysis](13-vibrational/) | `H2O_vib.inp` | Frequency analysis and IR spectrum |
| 14 | [TDDFT](14-tddft/) | `ethylene_tddft.inp` | Time-dependent DFT for excited states |

### ⭐⭐⭐⭐⭐ Expert

| # | Chapter | Input Files | Description |
|---|---------|-------------|-------------|
| 15 | [Hybrid Functionals](15-mp2/) | `water_mp2.inp` | PBE0 hybrid functional calculation |
| 16 | [Band Structure](16-gw-bse/) | `silicon_gw.inp` | Electronic band structure of silicon |

## Project Structure

```
cp2k_learn/
├── 01-hello-cp2k/          # Chapter 1: First CP2K run
├── 02-water-static/        # Chapter 2: Static DFT
├── 03-basis-sets/           # Chapter 3: Basis set comparison
├── 04-cutoff-convergence/  # Chapter 4: Cutoff convergence
├── 05-periodic-silicon/    # Chapter 5: Periodic systems
├── 06-geo-opt/             # Chapter 6: Geometry optimization
├── 07-cell-opt/            # Chapter 7: Cell optimization
├── 08-md-nve/              # Chapter 8: NVE molecular dynamics
├── 09-md-thermostats/      # Chapter 9: NVT/NPT thermostats
├── 10-metadynamics/        # Chapter 10: Metadynamics
├── 11-xtb/                 # Chapter 11: xTB semi-empirical
├── 12-qmmm-basics/        # Chapter 12: QM/MM hybrid
├── 13-vibrational/         # Chapter 13: Vibrational analysis
├── 14-tddft/               # Chapter 14: TDDFT excited states
├── 15-mp2/                 # Chapter 15: Hybrid functionals
├── 16-gw-bse/              # Chapter 16: Band structure
├── index.html              # Interactive web interface (SPA)
├── app.js                  # Web interface application logic
├── styles.css              # Web interface styling
├── LICENSE                 # GPL v3 License
└── README.md               # This file
```

Each chapter directory contains:
- `README.md` — Tutorial documentation (in Chinese)
- `*.inp` — CP2K input file(s)
- `*.xyz` — Molecular/crystal structure file(s) (where applicable)

## Web Interface

The project includes a single-page application for browsing the tutorials interactively. Open `index.html` in any modern browser.

Features:
- Chapter navigation with difficulty indicators
- Progress tracking (stored in localStorage)
- Syntax-highlighted input file viewer
- Responsive design for mobile and desktop

## Contributing

Contributions are welcome! Here's how you can help:

1. **Report bugs** — Open an issue describing the problem with the CP2K input or documentation
2. **Suggest new chapters** — Propose topics not yet covered (e.g., AIMD, NEB, RPA)
3. **Improve documentation** — Fix typos, add explanations, improve translations
4. **Add examples** — Submit new tested `.inp` files following the existing chapter structure

### Development Workflow

```bash
# Fork and clone
git fork https://github.com/user/cp2k-learn.git
git clone https://github.com/YOUR_USERNAME/cp2k-learn.git

# Create a feature branch
git checkout -b feature/new-chapter

# Make changes and test
cd 17-new-chapter
cp2k.ssmp -i new_calculation.inp

# Commit and push
git add .
git commit -m "Add Chapter 17: New Topic"
git push origin feature/new-chapter

# Open a Pull Request
```

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

---

*All CP2K input files have been validated against CP2K v2025.2 (`cp2k.ssmp`).*
