# ArDCA_PPI_sampler – Modeling Protein-Protein Sequences with promiscuos paraog interactions using ArDCA Evolution

> Building trustworthy datasets for interacting paralogs—especially **promiscuous** interactions—is hard.  
> This project creates **synthetic yet realistic** interaction datasets using **ArDCA** (coevolutionary autoregressive models) and **phylogenetic sampling**, then benchmarks paralog matching algorithms.

---

## Table of Contents
- [Motivation](#motivation)
- [Highlights](#highlights)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
  - [ArDCA: conditional partner generation](#ardca-conditional-partner-generation)
  - [Promiscuity modes (1→many / many→1 / many↔many)](#promiscuity-modes-1many--many1--manymany)
  - [Phylogenetic sampling](#phylogenetic-sampling)
- [Evaluation: paralog matching](#evaluation-paralog-matching)
- [Repository Structure](#repository-structure)
- [Figures](#figures)
- [Roadmap](#roadmap)
- [Citing](#citing)
- [License](#license)

---

## Motivation
Annotated datasets for **interacting paralogs** are scarce; for **promiscuous** interactions they are practically nonexistent. Coevolution-based generative models have repeatedly produced sequences that **fold and function** like natural ones. Here, we leverage them to **generate labeled interaction datasets** (one-to-one and promiscuous) and to **quantify** how matching methods behave under different signals—**coevolution**, **phylogeny**, and their combination.

---

## Highlights
- **ArDCA** for **fast, conditional** sampling of interacting partners.
- **Phylogenetic sampling** to reintroduce evolutionary relatedness missing from equilibrium-like sampling.
- Flexible **interaction graphs**: one-to-one, **one-to-many**, **many-to-one**, **many-to-many**.
- Reproducible **benchmarks** for matching algorithms (e.g., IPA) with TP/FP tracking and recruitment dynamics.

---

## Installation

### Requirements
- **Julia ≥ 1.10**
- `git`
- (Optional) Python + Matplotlib if exporting plots via PyPlot.

### Setup
```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>
julia --project -e 'using Pkg; Pkg.instantiate()'
