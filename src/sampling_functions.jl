#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
	sampling_one2one(inputfastafile::String, outputfastafile::String, M::Int64, La::Int64; kwds...) -> (Matrix{Int64}) 
Generates one-to-one samples using ArDCA models from natural one-to-one interacting sequences.
It returns the generated samples in 𝐿×𝑀 matrix of integer and saves them in a H5DF file.


* 'inputfastafile': is an input file in fasta format (for the natural one-to-one interacting sequences),
* 'outputfastafile': is an output file in fasta format (for generated one-to-one samples), 
* '𝑀': is the number of generated samples.

The keyword arguments for the sampling method with their default value:

* 'lambdaJ'=0.0001 coupling L₂ regularization parameter (Lagrange multiplier),
* 'lambdaH'=0.000001 field L₂ regularization parameter (Lagrange multiplier),
* 'permorder'="Natural" is the autoregressive order. 

Example of use:

inputfastafile = "./data/HKa-RRa_for_arDCA.fasta"
outputfastafile = "./data/HK-RR_ArDCA_Natural_1to1_M=10000.h5"

sampling_one2one(inputfastafile, outputfastafile, 10000, 64)
"""

function sampling_one2one(inputfastafile::String, 
		outputfastafile::String, 
		M::Int64, 
		La::Int64; 
		lambdaJ_s=10^(-4),
		lambdaH_s=10^(-6), 
		permorder_s=:NATURAL)

    #inputfastafile is an input file in fasta format.
    #outputfastafile is an output file in fasta format.
    #M is the number of generated samples.
    #La is the number of amino acids for family A.
    #lambdaJ=0.0001 coupling L₂ regularization parameter (Lagrange multiplier)
    #lambdaH=0.000001 field L₂ regularization parameter (Lagrange multiplier)
    #permorder is the autoregressive order. 

    #---------------------------------------------------------------------------------------
    # Learns the field and the coupling parameters h, J from the MSA using the ArDCA method...
    #...parameters are stored in arnet and the other algorithms variables are stored in arvar...
    #...lambdaJ=10^(-4),lambdaH=10^(-6), permorder
    arnet, arvar = ardca(inputfastafile, verbose=false, lambdaJ=lambdaJ_s, lambdaH=lambdaH_s, permorder=permorder_s)

    #---------------------------------------------------------------------------------------
    # Generates one-to-one 𝑀 sequences using the 'sample_with_weights' method from ArDCA package...
    #...'generated_alignment' is a 𝐿×𝑀 matrix of integer where 𝐿 is the alignment's length. 
    _, generated_alignment = sample_with_weights(arnet, M)
    interact = hcat(1:M,1:M) #interaction map for one-to-one interactions.

    #---------------------------------------------------------------------------------------
    # Writes generated sequences into 'outputfastafile' as a HDF5 file.  
    write_h5(generated_alignment[1:La,:], generated_alignment[La+1:end,:], outputfastafile, interact)

    return generated_alignment

end


#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
	sampling_one2many_known(inputfastafile::String, outputfastafile::String, Ma::Int64, La::Int64, λ::Int64; kwds...) -> (Matrix{Int64})
Generates one-to-many samples using ArDCA models from natural one-to-one interacting sequences. For each family A
sequence, there are λ interacting partners from family B. It returns the generated samples in 𝐿×𝑀b matrix of 
integer, where 𝑀b = λ×𝑀a and Ma is the number of sampled A sequences. Generated MSAs are saved in a H5DF file.


* 'inputfastafile': is the input file in fasta format (for the natural one-to-one interacting sequences),
* 'outputfastafile': is the output file in fasta format (for generated one-to-λ samples), 
* '𝑀a': is the number of unique samples for family A.
* 'La': is the number of amino acids for family A.
* 'λ': is the number of interacting partners for each A sequence.

The keyword arguments for the sampling method with their default value:

* 'lambdaJ'=0.0001 coupling L₂ regularization parameter (Lagrange multiplier),
* 'lambdaH'=0.000001 field L₂ regularization parameter (Lagrange multiplier),
* 'permorder'="Natural" is the autoregressive order. 

Example of use:

inputfastafile = "./data/HKa-RRa_for_arDCA.fasta"
outputfastafile = "./data/HK-RR_ArDCA_Natural_1to2_M=5000.h5"

sampling_one2many_known(inputfastafile, outputfastafile, 5000, 64, 2)
"""

function sampling_one2many_known(inputfastafile::String, 
		outputfastafile::String, 
		Ma::Int64, 
		La::Int64, 
		λ::Int64; 
		lambdaJ_s=10^(-4),
		lambdaH_s=10^(-6), 
		permorder_s=:NATURAL)

    #inputfastafile is an input file in fasta format.
    #outputfastafile is an output file in fasta format.
    #𝑀a is the number of unique samples for family A.
    #La is the number of amino acids for family A.
    #λ is the known number of interaction partners for each A sequence. 
    #lambdaJ=0.0001 coupling L₂ regularization parameter (Lagrange multiplier)
    #lambdaH=0.000001 field L₂ regularization parameter (Lagrange multiplier)
    #permorder is the autoregressive order. 

    #---------------------------------------------------------------------------------------
    # Learns the field and the coupling parameters h, J from the MSA using the ArDCA method...
    #...parameters are stored in arnet, and the other algorithms variables are stored in arvar...
    #...lambdaJ=10^(-4),lambdaH=10^(-6), permorder
    arnet, arvar = ardca(inputfastafile, verbose=false, lambdaJ=lambdaJ_s, lambdaH=lambdaH_s, permorder=permorder_s)

    #---------------------------------------------------------------------------------------
    # Generates one-to-one 𝑀 sequences using the 'sample_with_weights' method from ArDCA package...
    #...'generated_alignment' is a 𝐿×𝑀 matrix of integer where 𝐿 is the alignment's length. 
    _, generated_alignment_1to1 = sample_with_weights(arnet, Ma)

    #---------------------------------------------------------------------------------------
    # Generates one-to-λ
    L = size(generated_alignment_1to1, 1)
    generated_alignment_1toλ = fill(0, L, λ * Ma)
    interact = Array{Int64,2}(undef, λ * Ma, 2)
    
    for i in 1:Ma
        _, generated_alignment_1toλ[:, (i-1)*λ + 1:i*λ] = sample_subsequence(generated_alignment_1to1[1:La,i], arnet, λ)
        # We save interaction map for one-to-λ interactions.
        for j in (i-1)*λ + 1:i*λ
            interact[j, 1] = i
            interact[j, 2] = j
        end
    end

    #---------------------------------------------------------------------------------------
    # Writes generated sequences into 'outputfastafile' as a HDF5 file.  
    write_h5(generated_alignment_1to1[1:La,:], generated_alignment_1toλ[La+1:end,:], outputfastafile, interact)

    return generated_alignment_1toλ

end


#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
	sampling_one2many_unknown(inputfastafile::String, outputfastafile::String, Ma::Int64, La::Int64, λ::Int64; kwds...) -> (Matrix{Int64})
Generates one-to-many samples using ArDCA models from natural one-to-one interacting sequences. For each family A 
sequence, there are λi interacting partners from family B, λi is sampled from a Poisson distribution with mean λ. 
It returns the generated samples in 𝐿×𝑀 matrix of integer, where 𝑀 is the number of sampled sequences. Generated 
MSAs are saved in a H5DF file:

* 'inputfastafile': is the input file in fasta format (for the natural one-to-one interacting sequences),
* 'outputfastafile': is the output file in fasta format (for generated one-to-λ samples), 
* '𝑀a': is the number of unique samples for family A.
* 'La': is the number of amino acids for family A. 
* 'λ': is the mean of Poisson distribution (used for sampling interaction partners for each A sequence).

The keyword arguments for the sampling method with their default value:

* 'lambdaJ'=0.0001 coupling L₂ regularization parameter (Lagrange multiplier),
* 'lambdaH'=0.000001 field L₂ regularization parameter (Lagrange multiplier),
* 'permorder'="Natural" is the autoregressive order. 

Example of use:

inputfastafile = "./data/HKa-RRa_for_arDCA.fasta"
outputfastafile = "./data/HK-RR_ArDCA_Natural_one2many_Posisson_mu=2_M=5000.h5"

sampling_one2many_known(inputfastafile, outputfastafile, 5000, 64, 2)
"""

function sampling_one2many_unknown(inputfastafile::String, 
		outputfastafile::String, 
		Ma::Int64, 
		La::Int64, 
		λ::Int64; 
		lambdaJ_s=10^(-4), 
		lambdaH_s=10^(-6), 
		permorder_s=:NATURAL)

    #inputfastafile is an input file in fasta format.
    #outputfastafile is an output file in fasta format.
    #𝑀a is the number of generated samples.
    #La is the number of amino acids for family A.
    #λ is the mean of Poisson distribution (used for sampling interaction partners for each A sequence).
    #lambdaJ=0.0001 coupling L₂ regularization parameter (Lagrange multiplier).
    #lambdaH=0.000001 field L₂ regularization parameter (Lagrange multiplier).
    #permorder is the autoregressive order. 

    #---------------------------------------------------------------------------------------
    # Learns the field and the coupling parameters h, J from the MSA using the ArDCA method...
    #...parameters are stored in arnet, and the other algorithms variables are stored in arvar...
    #...lambdaJ=10^(-4),lambdaH=10^(-6), permorder
    arnet, arvar = ardca(inputfastafile, verbose=false, lambdaJ=lambdaJ_s, lambdaH=lambdaH_s, permorder=permorder_s)

    #---------------------------------------------------------------------------------------
    # Generates one-to-one 𝑀 sequences using the 'sample_with_weights' method from ArDCA package...
    #...'generated_alignment' is a 𝐿×𝑀 matrix of integer where 𝐿 is the alignment's length. 
    _, generated_alignment_1to1 = sample_with_weights(arnet, Ma)

    #---------------------------------------------------------------------------------------
    # Samples interaction partner for each A sequence from a Poisson distribution with mean λ
    inter_partners = rand(Poisson(λ), Ma)  # Generate Poisson samples
    inter_partners = max.(inter_partners, 1)  # Replace zeros with ones
    Mb = sum(inter_partners) 
    
    #---------------------------------------------------------------------------------------
    # Generates one-to-many following 'inter_partners' (the Poisson distribution with mean λ)
    L = size(generated_alignment_1to1, 1)
    generated_alignment_1tomany = fill(0, L, Mb)
    interact = Array{Int64, 2}(undef, Mb, 2)
    count_seq = 0
    
    for i in 1:Ma
        λi = inter_partners[i]
        _, generated_alignment_1tomany[:, count_seq + 1:count_seq + λi] = sample_subsequence(generated_alignment_1to1[1:La,i], arnet, λi)
        # We save interaction map for one-to-many interactions.
        for j in count_seq + 1:count_seq + λi
            interact[j, 1] = i
            interact[j, 2] = j
        end
        count_seq += λi
    end

    #---------------------------------------------------------------------------------------
    # Writes generated sequences into 'outputfastafile' as a HDF5 file.  
    write_h5(generated_alignment_1to1[1:La,:], generated_alignment_1tomany[La+1:end,:], outputfastafile, interact)

    return generated_alignment_1tomany

end



#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
	sampling_many2one_known(inputfastafile::String, outputfastafile::String, Mb::Int64, Lb::Int64, λ::Int64; kwds...) -> (Matrix{Int64})
Generates many-to-one samples using ArDCA models from natural one-to-one interacting sequences. For each family 
B sequence, there are λ interacting partners from family A. It returns the generated samples in 𝐿×𝑀a matrix of 
integer, where 𝑀a = λ×𝑀b and Mb is the number of sampled B sequences. Generated MSAs are saved in a H5DF file:

* 'inputfastafile': is the input file in fasta format (for the natural one-to-one interacting sequences),
* 'outputfastafile': is the output file in fasta format (for generated λ-to-one samples), 
* '𝑀b': is the number of unique samples for family B.
* 'Lb': is the number of amino acids for family B.
* 'λ': is the number of interacting partners for each B sequence.

The keyword arguments for the sampling method with their default value:

* 'lambdaJ'=0.0001 coupling L₂ regularization parameter (Lagrange multiplier),
* 'lambdaH'=0.000001 field L₂ regularization parameter (Lagrange multiplier),
* 'permorder'="Natural" is the autoregressive order. 

Example of use:

inputfastafile = "./data/HKa-RRa_for_arDCA.fasta"
outputfastafile = "./data/HK-RR_ArDCA_antiNatural_2to1_M=5000.h5"

sampling_many2one_known(inputfastafile, outputfastafile, 5000, 112, 2)
"""

function sampling_many2one_known(inputfastafile::String, 
		outputfastafile::String, 
		Mb::Int64, 
		Lb::Int64, 
		λ::Int64; 
		lambdaJ_s=10^(-4), 
		lambdaH_s=10^(-6), 
		permorder_s=sort(1:176, rev=true))

    #inputfastafile is an input file in fasta format.
    #outputfastafile is an output file in fasta format.
    #𝑀b is the number of unique samples for family B.
    #Lb is the number of amino acids for family B.
    #λ known number of interaction partners for each B sequence.
    #lambdaJ=0.0001 coupling L₂ regularization parameter (Lagrange multiplier)
    #lambdaH=0.000001 field L₂ regularization parameter (Lagrange multiplier)
    #permorder is the autoregressive order. 

    #---------------------------------------------------------------------------------------
    # Learns the field and the coupling parameters h, J from the MSA using the ArDCA method...
    #...parameters are stored in arnet, and the other algorithms variables are stored in arvar...
    #...lambdaJ=10^(-4),lambdaH=10^(-6), permorder
    arnet, arvar = ardca(inputfastafile, verbose=false, lambdaJ=lambdaJ_s,lambdaH=lambdaH_s, permorder=permorder_s)

    #---------------------------------------------------------------------------------------
    # Generates one-to-one 𝑀 sequences using the 'sample_with_weights' method from ArDCA package...
    #...'generated_alignment' is a 𝐿×𝑀 matrix of integer where 𝐿 is the alignment's length. 
    _, generated_alignment_1to1 = sample_with_weights(arnet, Mb)

    #---------------------------------------------------------------------------------------
    # Generates many-to-one
    L = size(generated_alignment_1to1, 1)
    generated_alignment_λto1 = fill(0, L, λ * Mb)
    interact = Array{Int64, 2}(undef, λ * Mb, 2)
    
    for i in 1:Mb
        generated_alignment_λto1[:, (i-1)*λ + 1:i*λ] = sample_subsequence_many2one(generated_alignment_1to1[L-Lb+1:L,i], arnet, λ)
        # We save interaction map for λ-to-one interactions.
        for j in (i-1)*λ + 1:i*λ
            interact[j, 1] = j
            interact[j, 2] = i
        end
    end

    #---------------------------------------------------------------------------------------
    # Writes generated sequences into 'outputfastafile' as a HDF5 file.
    La = L - Lb
    println("L = ", L, ", La = ", La, ", Lb = ", Lb)
    write_h5(generated_alignment_λto1[1:La,:], generated_alignment_1to1[La+1:L,:], outputfastafile, interact)

    return generated_alignment_λto1

end



#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
    sampling_many2one_unknown(inputfastafile::String, outputfastafile::String, Mb::Int64, Lb::Int64, λ::Int64; kwds...) -> (Matrix{Int64})
Generates many-to-one samples using ArDCA models from natural one-to-one interacting sequences.
For each family B sequence, there are λi interacting partners from family A, λi is sampled from
a Poisson distribution with mean λ. It returns the generated samples in 𝐿×𝑀 matrix of integer, 
where 𝑀 is the number of sampled sequences. Generated MSAs are saved in a H5DF file.

* 'inputfastafile': is the input file in fasta format (for the natural one-to-one interacting sequences),
* 'outputfastafile': is the output file in fasta format (for generated one-to-λ samples), 
* '𝑀b': is the number of unique samples for family B.
* 'Lb': is the number of amino acids for family B.
* 'λ': is the mean of Poisson distribution (used for sampling interaction partners for each A sequence).

The keyword arguments for the sampling method with their default value:

* 'lambdaJ'=0.0001 coupling L₂ regularization parameter (Lagrange multiplier),
* 'lambdaH'=0.000001 field L₂ regularization parameter (Lagrange multiplier),
* 'permorder'="Natural" is the autoregressive order. 

Example of use:

inputfastafile = "./data/HKa-RRa_for_arDCA.fasta"
outputfastafile = "./data/HHK-RR_ArDCA_antiNatural_many2one_Posisson_mu=2_M=5000.h5"

sampling_many2one_unknown(inputfastafile, outputfastafile, 5000, 112, 2)
"""


function sampling_many2one_unknown(inputfastafile::String, 
		outputfastafile::String, 
		Mb::Int64, 
		Lb::Int64, 
		λ::Int64; 
		lambdaJ_s=10^(-4), 
		lambdaH_s=10^(-6), 
		permorder_s=sort(1:176, rev=true))

    #inputfastafile is an input file in fasta format.
    #outputfastafile is an output file in fasta format.
    #𝑀b is the number of unique samples for family B.
    #Lb is the number of amino acids for family B.
    #λ is the mean of Poisson distribution (used for sampling interaction partners for each A sequence).
    #lambdaJ=0.0001 coupling L₂ regularization parameter (Lagrange multiplier)
    #lambdaH=0.000001 field L₂ regularization parameter (Lagrange multiplier)
    #permorder is the autoregressive order. 

    #---------------------------------------------------------------------------------------
    # Learns the field and the coupling parameters h, J from the MSA using the ArDCA method...
    #...parameters are stored in arnet, and the other algorithms variables are stored in arvar...
    arnet, arvar = ardca(inputfastafile, verbose=false, lambdaJ=lambdaJ_s,lambdaH=lambdaH_s, permorder=permorder_s)

    #---------------------------------------------------------------------------------------
    # Generate one-to-one 𝑀 sequences using the 'sample_with_weights' method from ArDCA package...
    #...'generated_alignment' is a 𝐿×𝑀 matrix of integer where 𝐿 is the alignment's length. 
    _, generated_alignment_1to1 = sample_with_weights(arnet, Mb)
    
    #---------------------------------------------------------------------------------------
    # Samples interaction partner for each A sequence from a Poisson distribution with mean λ
    inter_partners = rand(Poisson(λ), Mb)  # Generate Poisson samples
    inter_partners = max.(inter_partners, 1)  # Replace zeros with ones
    Ma = sum(inter_partners) 
    
    #---------------------------------------------------------------------------------------
    # Generates many-to-one following 'inter_partners' (the Poisson distribution with mean λ)
    L = size(generated_alignment_1to1, 1)
    generated_alignment_manyto1 = fill(0, L, Ma)
    interact = Array{Int64, 2}(undef, Ma, 2)
    count_seq = 0
    
    for i in 1:Mb
        λi = inter_partners[i]
        generated_alignment_manyto1[:, count_seq + 1:count_seq + λi] = sample_subsequence_many2one(generated_alignment_1to1[L-Lb+1:L,i], arnet, λi)
        # We save interaction map for many-to-one interactions.
        for j in count_seq + 1:count_seq + λi
            interact[j, 1] = j
            interact[j, 2] = i
        end
        count_seq += λi
    end

    #---------------------------------------------------------------------------------------
    # Writes generated sequences into 'outputfastafile' as a HDF5 file.
    La = L - Lb
    println(La)
    write_h5(generated_alignment_manyto1[1:La,:], generated_alignment_1to1[La+1:end,:], outputfastafile, interact)

    return generated_alignment_manyto1

end


#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
    sampling_many2many_unknown(inputfastafile::String, outputfastafile::String, Ma::Int64, La::Int64, λ::Int64; kwds...) -> (Matrix{Int64})
Generates many-to-many samples using ArDCA models from natural one-to-one interacting sequences.
For each family A sequence, there are λ interacting partners from family B. Then, from these λ 
related B sequences, we sample λ sequences from family A. It returns the generated samples in 
𝐿×𝑀 matrix of integer, where 𝑀 is the number of sampled sequences. Generated MSAs are saved in 
a H5DF file.


* 'inputfastafile': is the input file in fasta format (for the natural one-to-one interacting sequences),
* 'outputfastafile': is the output file in fasta format (for generated one-to-λ samples), 
* '𝑀b': is the number of unique samples for family B.
* 'Lb': is the number of amino acids for family B.
* 'λ': is the mean of Poisson distribution (used for sampling interaction partners for each A sequence).

The keyword arguments for the sampling method with their default value:

* 'lambdaJ'=0.0001 coupling L₂ regularization parameter (Lagrange multiplier),
* 'lambdaH'=0.000001 field L₂ regularization parameter (Lagrange multiplier).

Example of use:

inputfastafile = "./data/HKa-RRa_for_arDCA.fasta"
outputfastafile = "./data/HHK-RR_ArDCA_2to2_M=10000.h5"

sampling_many2many_known(inputfastafile, outputfastafile, 10000, 64, 2)
"""


function sampling_many2many_known(inputfastafile::String, 
		outputfastafile::String, 
		Ma::Int64, 
		La::Int64, 
		λ::Int64; 
		lambdaJ_s=10^(-4), 
		lambdaH_s=10^(-6))

    #inputfastafile is an input file in fasta format.
    #outputfastafile is an output file in fasta format.
    #𝑀a is the number of unique samples for family B.
    #La is the number of amino acids for family B.
    #λ is the mean of Poisson distribution (used for sampling interaction partners for each A sequence).
    #lambdaJ=0.0001 coupling L₂ regularization parameter (Lagrange multiplier)
    #lambdaH=0.000001 field L₂ regularization parameter (Lagrange multiplier) 

    #---------------------------------------------------------------------------------------
    # Learns the field and the coupling parameters h, J from the MSA using the ArDCA method...
    #...parameters are stored in arnet, and the other algorithms variables are stored in arvar...
    #...We use natural order
    arnet, arvar = ardca(inputfastafile, verbose=false, lambdaJ=lambdaJ_s, lambdaH=lambdaH_s, permorder=:NATURAL)

    #---------------------------------------------------------------------------------------
    # Generates one-to-one 𝑀 sequences using the 'sample_with_weights' method from ArDCA package...
    #...'generated_alignment' is a 𝐿×𝑀 matrix of integer where 𝐿 is the alignment's length. 
    _, generated_alignment_1to1 = sample_with_weights(arnet, Ma)

    #---------------------------------------------------------------------------------------
    # Generate one-to-λ
    L = size(generated_alignment_1to1, 1)
    generated_alignment_1toλ = fill(0, L, λ * Ma)
        
    for i in 1:Ma
        _, generated_alignment_1toλ[:, (i-1)*λ + 1:i*λ] = sample_subsequence(generated_alignment_1to1[1:La,i], arnet, λ)
    end

    #---------------------------------------------------------------------------------------
    # Learns the field and the coupling parameters h, J from the MSA using the ArDCA method...
    #...parameters are stored in arnet, and the other algorithms variables are stored in arvar...
    #...We use reverse natural order
    arnet_r, arvar_r = ardca(inputfastafile, verbose=false, lambdaJ=lambdaJ_s, lambdaH=lambdaH_s, permorder=sort(1:176, rev=true))

    #---------------------------------------------------------------------------------------
    # Generates λ-to-λ
    generated_alignment_λtoλ = fill(0, La, λ * Ma)
    
    for i in 1:Ma
        generated_alignment_λtoλ[:, (i-1)*λ + 1:i*λ] = sample_subsequence_many2many(generated_alignment_1toλ[La+1:L,(i-1)*λ + 1:i*λ], arnet_r, λ)
    end 

    #---------------------------------------------------------------------------------------
    # Generates λ-to-λ known interactions
    interact = Array{Int64,2}(undef, (λ^2) * Ma, 2)
        
    for i = 1:λ * Ma
        m = ceil(Int, i/λ)
        for j = 1:λ
            interact[λ*(i-1)+j, 1] = i
            interact[λ*(i-1)+j, 2] = λ*(m-1)+j
        end
    end

    #---------------------------------------------------------------------------------------
    # Writes generated sequences into 'outputfastafile' as a HDF5 file.  
    write_h5(generated_alignment_λtoλ, generated_alignment_1toλ[La+1:L,:], outputfastafile, interact)

end




#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
