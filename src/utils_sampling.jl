#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
    'sample_subsequence_many2one'
Returns a generated co-alignment in the form of a `L × msamples`. The co-alignment is forced to 
end with a sequence `x0` and then autoregressively generated:

* 'x0': is a B sequence fixed to generated 'msamples' A sequences, 
* 'arnet': is an ArDCA model learned from natural AB sequences in anti-Natural order,
* 'msamples': is the number of A samples generated from 'x0'.
"""


function sample_subsequence_many2one(x0::Vector{T}, 
		arnet::ArNet, 
		msamples::Int64) where {T<:Integer}

	#x0 is a B sequence fixed to generated 'msamples' A sequences, 
	#arnet is an ArDCA model learned from natural AB sequences in anti-Natural order,
	#msamples is the number of A samples generated from 'x0'.
    
    @extract arnet:H J p0 idxperm
    N = length(idxperm)
    
    #---------------------------------------------------------------------------------------
    # Checks sequences' length
    length(x0) < N || error("Subsequence too long for the model")    
    # Check right one-hot encoding
    all(x -> 1 ≤ x ≤ 21, x0) || error("Subsequence numeric code should be in 1..21 ")
    # Checks 'idxperm' is in anti-natural order.
    
    #---------------------------------------------------------------------------------------
    # Defines parameters
    l0 = length(x0)
    q = length(p0) 
    N = length(H) # here L is L-1 !!
    backorder = sortperm(idxperm);

    #---------------------------------------------------------------------------------------
    # Generates 'msamples' 
    res = Matrix{Int}(undef, N + 1, msamples)

    Threads.@threads for i in 1:msamples
        totH = Vector{Float64}(undef, q)
        sample_z = -ones(Int, N + 1)
        @turbo for k in 1:l0
            sample_z[k] = x0[l0 + 1 - k]
        end

        for site in 1:N
            Js = J[site]
            h = H[site]
            copy!(totH, h)
            @turbo for i in 1:site
                for a in 1:q
                    totH[a] += Js[a, sample_z[i], i]
                end
            end
            p = softmax(totH)
            if sample_z[site+1] == -1
                sample_z[site+1] = wsample(1:q, p)
            end
        end
        res[:, i] .= sample_z        
    end
    
    return permuterow!(res, backorder)
    
end



#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
"""
	'sample_subsequence_many2many'
Returns a generated alignment in the form of a `La × msamples`. The alignment is forced to
be related with a sequence `x0m` and then autoregressively generated.

* 'x0m': is a matrix of B sequences fixed to generated 'msamples' of A sequences, 
* 'arnet': is an ArDCA model learned from natural AB sequences in anti-Natural order,
* 'msamples': is the number of A samples generated from 'x0m'.
"""


function sample_subsequence_many2many(x0m::Matrix{T}, 
		arnet::ArNet, 
		msamples::Int64) where {T<:Integer}
	
	#x0m is a matrix of B sequences fixed to generated 'msamples' of A sequences, 
	#arnet is an ArDCA model learned from natural AB sequences in anti-Natural order,
	#msamples is the number of A samples generated from 'x0m'.
    
    @extract arnet:H J p0 idxperm
    N = length(idxperm)
    
    #---------------------------------------------------------------------------------------
    # Checks sequences' length
    l0, lambda = size(x0m)
    l0 < N || error("Subsequence too long for the model")
    # Checks right one-hot encoding
    all(x -> 1 ≤ x ≤ 21, x0m) || error("Subsequences numeric code should be in 1..21 ")
    # Checks 'idxperm' is in anti-natural order.
    
    #---------------------------------------------------------------------------------------
    # Defines parameters
    q = length(p0) 
    N = length(H) # here N is N-1 !!
    backorder = sortperm(idxperm)

    #---------------------------------------------------------------------------------------
    # Generates 'msamples' 
    Lsamples = N - l0 + 1
    res = Matrix{Int}(undef, Lsamples, msamples)

    Threads.@threads for i in 1:msamples
        totH = Vector{Float64}(undef, q)        
        sample_z = -ones(Int, N + 1, lambda)
        @turbo for l in 1:lambda, k in 1:l0
            sample_z[k,l] = x0m[l0 + 1 - k, l]
        end

        for site in 1:N
            Js = J[site]
            h = H[site]
            copy!(totH, h)
            @turbo for l in 1:lambda
                for i in 1:site
                    for a in 1:q
                        totH[a] += Js[a, sample_z[i,l], i]/lambda
                    end
                end
            end
            p = softmax(totH)
            if sample_z[site+1,1] == -1
                sample_aa = wsample(1:q, p)
                @turbo for l in 1:lambda
                    sample_z[site+1,l] = sample_aa
                end
            end
        end
        res[:,i] = sample_z[l0+1:N+1, 1]
    end
    
    return permuterow!(res, backorder[l0+1:N+1])
    
end



#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------
