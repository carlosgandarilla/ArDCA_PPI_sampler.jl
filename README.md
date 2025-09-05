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
julia --project -e 'using Pkg; Pkg.instantiate()
```

---

## Quick Start

### 1) Train/load ArDCA and generate one-to-one pairs
```julia
using ProjectName  # replace with your module name

m = fit_ardca("data/msa_A.fasta", "data/msa_B.fasta"; order=:natural)  # or :entropic
pairs = generate_pairs(m; N=4000, out="results/pairs_one_to_one.fasta")
```

### 2) Generate one-to-many with Poisson multiplicity
```julia
ds = generate_one_to_many(m; N_A=2000, lambda=2.0, out="results/one_to_many/")
```

### 3) Generate many-to-many on a specified bipartite graph
```julia
g  = load_graph("configs/graph.toml")   # degrees, connected components
mm = generate_many_to_many(m, g; out="results/many_to_many/")
```

---

## How It Works

### ArDCA: conditional partner generation
ArDCA models a protein sequence $\mathbf{a}=(a_1,\dots,a_L)$ with an autoregressive factorization:

$$
P(\mathbf{a})=\prod\limits_{i=1}^{L} P\big(a_i \mid a_{i-1},\ldots, a_1\big).
$$

Each conditional is parameterized in softmax form:

$$
P\big(a_i \mid a_{i-1},\ldots, a_1\big)=
\frac{\exp \left \{ h_i(a_i) + \sum\limits_{j=1}^{i-1} J_{ij}\big(a_i,a_j\big) \right\} }{Z_i\big(a_1,\ldots,a_{i-1}\big)}.
$$

For interacting families $A$ and $B$, we generate partners by conditioning:

$$
P(\mathbf{a},\mathbf{b})=P(\mathbf{a})\,P(\mathbf{b}\mid \mathbf{a}),\qquad
P(\mathbf{b}\mid \mathbf{a})=\prod_{i=1}^{L_B} P\!\big(b_i \mid b_1,\ldots,b_{i-1},\mathbf{a}\big),
$$

with cross-family couplings:

$$
P\!\big(b_i \mid b_1,\ldots,b_{i-1},\mathbf{a}\big)=
\frac{\exp\!\left\{h_i(b_i)+\sum_{j=1}^{i-1}J_{ij}\!\big(b_i,b_j\big)+\sum_{k=1}^{L_A} L_{ik}\!\big(b_i,a_k\big)\right\}}
     {Z_i\!\big(b_1,\ldots,b_{i-1},\mathbf{a}\big)}.
$$

**Notes**
- **Positional order**: “natural” vs “entropic.” Entropic often fits better; natural is convenient for pairing workflows.
- **Sampling temperature**: optionally sample at inverse temperature $\beta>1$ to bias toward lower energies (default $\beta=1$).

### Promiscuity modes (1→many / many→1 / many↔many)
- **One-to-many**: for each $\mathbf{a}^m$, sample $c_m$ partners in $B$.
  - Known multiplicity: fixed $c_m$.
  - Unknown multiplicity: draw $c_m$ from

    $$
    P(c_m)=\frac{\lambda^{c_m}e^{-\lambda}}{c_m!}.
    $$

- **Many-to-one**: generate one $\mathbf{b}$ conditioned on multiple $\{\mathbf{a}_\alpha\}$.
  - Either average encodings of the $\mathbf{a}_\alpha$ (simpler), or concatenate them (stronger, may need per-$M$ models).
- **Many-to-many**: define a bipartite interaction graph (node degrees), then alternate conditional sampling until the component is filled.

### Phylogenetic sampling
ArDCA alone lacks phylogenetic structure. To reintroduce realistic relatedness:
1. Pick a **root** natural sequence.  
2. For each branch (length $t$), propose $L\,\mu\,t$ mutations.  
3. Accept/reject using a Potts-model energy (parameters learned on the natural co-MSA via BM-DCA).  
4. Traverse the tree to obtain sequences mirroring phylogenetic relationships.

**Practical tip (HK–RR case)**: $\mu \approx 1.8$ best matches the low-distance tail of natural Hamming-distance distributions.

---

## Evaluation: paralog matching
We benchmark with **IPA** (Iterative Paralog Matching), tracking:
- True-positive / false-positive fractions across iterations
- Recruitment dynamics (how new pairs enter the concatenated alignment)
- Sensitivity to **paralog multiplicity** and **dataset size**

Expected trends:
- **Coevolution without phylogeny** ⇒ highest matching accuracy  
- **Phylogeny alone** ⇒ can reinforce false positives  
- **Natural datasets** (both signals) ⇒ intermediate performance

**Example**
```julia
r = run_ipa("results/pairs_one_to_one.fasta";
            nincrement=6, n_runs=50, split_by_species=true)

plots = plot_recruitment_dynamics(r; out="results/plots/")
```

---

## Repository Structure
```
.
├─ src/                      # Julia code: ArDCA, sampling, graphs, IPA wrappers, plotting
├─ data/                     # MSAs, phylo trees, root sequences (placeholders)
├─ configs/                  # Graph specs (TOML/JSON)
├─ results/                  # Generated datasets and plots
├─ docs/
│  └─ figures/               # Images used in README
└─ README.md
```

---

## Figures
- `docs/figures/ardca_schematic.png` — ArDCA conditional structure  
- `docs/figures/generation_modes.png` — one→many / many→1 / many↔many flows  
- `docs/figures/phylo_hist.png` — Hamming-distance distributions vs $\mu$

> Replace these paths with your actual file names (or leave as TODOs).

---

## Roadmap
- Graph-conditional training to avoid per-$M$ retraining in many-to-one
- Larger families (including eukaryotic) and richer promiscuity patterns
- Side-by-side: IPA vs MI vs transformer-based pairing on these datasets

---

## Citing
If you use this repository or datasets, please cite:

```
@misc{paralog_generation_2025,
  title  = {Interplay of Phylogeny and Coevolution in Partner Selection:
            Labeled datasets for paralog matching with ArDCA and phylogenetic sampling},
  author = {<Your Name> and <Collaborators>},
  year   = {2025},
  note   = {GitHub repository},
  url    = {https://github.com/<user>/<repo>}
}
```

Also consider citing ArDCA and related coevolutionary work as appropriate.

---

## License
Add a LICENSE file (e.g., MIT or Apache-2.0) and reference it here.

