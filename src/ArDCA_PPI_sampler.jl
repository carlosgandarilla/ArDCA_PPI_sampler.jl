#begin of module
module ArDCA_PPI_sampler

	#--------------------------------------------------------------------------------------
	#import	
	
	#import ArDCA: ardca, ArNet, sample_with_weights, sample_subsequence, permuterow!, wsample
	import ExtractMacro: @extract
	import DCAUtils: read_fasta_alignment, compute_weighted_frequencies, compute_weights, add_pseudocount
	import LoopVectorization: @turbo
	import HDF5: h5write, h5read
	import CurveFit: linear_fit
	import PyPlot: figure, hist, plot, scatter, xlim, ylim, xticks, yticks, xlabel, ylabel, xscale, yscale, legend, savefig
	import LaTeXStrings: @L_str
	import Distances: pairwise, Hamming
	import StatsBase: cor
	import Distributions: Poisson#, cor
	
	#--------------------------------------------------------------------------------------
	#export

	export sampling_one2one,
	       sampling_one2many_known,
	       sampling_one2many_unknown,
	       sampling_many2one_known,
	       sampling_many2one_unknown,
	       sampling_many2many_known,
	       computeCijk,
	       dij_hist_generative,
	       Cijk_fig_generative,
	       Cij_fig_generative,
	       fi_fig_generative
	       
	#--------------------------------------------------------------------------------------
	#include julia source files
	
	include("three_point_correlations.jl")
	include("sampling_generative_properties.jl")
	include("sampling_functions.jl")
	include("utils_sampling.jl")
	include("utils_data.jl")

end # end module
