```@meta
CurrentModule = ArDCA_PPI_sampler
```

# ArDCA_PPI_sampler.jl documentation

Understanding protein-protein interactions (PPI) is essential in bioinformatics, as they play
a key role in numerous biological processes. However, experimentally identifying these interactions
is challenging, particularly in promiscuous binding cases where proteins interact with multiple partners. 
For promiscuous interactions, no comprehensive annotated datasets are available, and for interacting 
paralogs are scarce. We have adapted ArDCA sequence generation techniques to simulate interacting 
sequences based on various interaction graphs. Because of its autoregressive nature, the model allows 
to sample sequences conditioned on each other. These labeled datasets can then be used to test 
paralog matching methods and evaluate their performance. 

We will assume we have a concatenated Multiple Sequence Alignment (co-MSA) in FASTA format which is 
obtained by concatenating known interaction partners from two protein families. We model the interaction
between proteins by statistically conditioning interacting sequences on each other using a 
coevolution-based method, i.e. by using pairwise models to sample sequences from joint-probability 
distributions learned on paired MSAs.

The accuracy of generative models can be assessed by comparing statistical properties of the sampled
sequences with those of natural sequences. Specifically, one-point frequencies, two-point correlations, 
and three-point correlations provide key metrics to evaluate how well the model captures the underlying 
sequence distributions. Additionally, the Hamming distance between sampled sequences can be analyzed to 
assess sequence diversity and similarity patterns. By examining these measures, we can determine the extent 
to which the generated sequences faithfully reproduce the characteristics of natural protein sequences.

We will use the histidine kinase (HK) and response regulators (RR) protein families available in `data/` 
folder of the package/

The code is written in [Julia](http://julialang.org). It requires Julia `>=1.10` or later.


## Installation

To install the package, use Julia's package manager: from the Julia REPL, type `]` to enter the Pkg
REPL mode and run:

```
(v1.?) pkg> add https://github.com/carlosgandarilla/ArDCA_PPI_sampler
```

Then load it with:

```
julia> using ArDCA_PPI_sampler
```

## Multithreading

To take full advantage of multicore parallel computing, Julia should be launched with

```
$ julia -t nthreads # put here nthreads equal to the number of cores you want to use
```

More information [here](https://docs.julialang.org/en/v1.6/manual/multi-threading/)

## Reference

```@docs
sampling_one2one
```
```@docs
sampling_one2many_known
```
```@docs
sampling_one2many_unknown
```
```@docs
sampling_many2one_known
```
```@docs
sampling_many2one_unknown
```
```@docs
sampling_many2many_known
```
```@docs
computeCijk
```
```@docs
dij_hist_generative
```
```@docs
Cijk_fig_generative
```
```@docs
Cij_fig_generative
```
```@docs
fi_fig_generative
```
