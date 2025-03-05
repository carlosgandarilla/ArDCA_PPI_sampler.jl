using ArDCA_PPI_sampler
using Test

@testset "MyJuliaPackage Tests" begin
    @test true  # Placeholder test to ensure test suite runs
end

@testset "Testing ArDCA_PPI_sampler" begin

    # Define input data
    inputfastafile = "./data/HKa-RRa_for_arDCA.fasta"
    
    # Test one-to-one sampling (check if it runs without errors)
    try
        #define path for generated one-to-one samples
        outputfastafile = "./data/HK-RR_ArDCA_Natural_1to1_M=10000.h5"
        sampling_one2one(inputfastafile, outputfastafile, 10000, 64)
        @test true   # If no error occurs, this test passes
    catch e
        @test false  # If an error occurs, this test fails
    end

    # Test one-to-many sampling (when we know how many)
    try
        #define path for generated one-to-two samples
        outputfastafile = "./data/HK-RR_ArDCA_Natural_1to2_M=5000.h5"
        sampling_one2many_known(inputfastafile, outputfastafile, 5000, 64, 2)
        @test true   # If no error occurs, this test passes
    catch e
        @test false  # If an error occurs, this test fails
    end

    # Test one-to-many sampling (when we do not know how many)
    try
        #define path for generated one-to-two samples
        outputfastafile = "./data/HK-RR_ArDCA_Natural_one2many_Posisson_mu=2_M=5000.h5"
        sampling_one2many_known(inputfastafile, outputfastafile, 5000, 64, 2)
        @test true   # If no error occurs, this test passes
    catch e
        @test false  # If an error occurs, this test fails
    end

    # Test many-to-one sampling (when we know how many)
    try
        #define path for generated one-to-two samples
        outputfastafile = "./data/HK-RR_ArDCA_antiNatural_2to1_M=5000.h5"
        sampling_many2one_known(inputfastafile, outputfastafile, 5000, 112, 2)
        @test true   # If no error occurs, this test passes
    catch e
        @test false  # If an error occurs, this test fails
    end

    # Test many-to-one sampling (when we do not know how many)
    try
        #define path for generated one-to-two samples
        outputfastafile = "./data/HHK-RR_ArDCA_antiNatural_many2one_Posisson_mu=2_M=5000.h5"
        sampling_many2one_unknown(inputfastafile, outputfastafile, 5000, 112, 2)
        @test true   # If no error occurs, this test passes
    catch e
        @test false  # If an error occurs, this test fails
    end

    # Test many-to-many sampling (when we know how many)
    try
        #define path for generated one-to-two samples
        outputfastafile = "./data/HHK-RR_ArDCA_2to2_M=10000.h5"
        sampling_many2many_known(inputfastafile, outputfastafile, 10000, 64, 2)
        @test true   # If no error occurs, this test passes
    catch e
        @test false  # If an error occurs, this test fails
    end







    

    # Test for save_hdf5 (check if it runs without errors)
    try
        save_hdf5([1, 2, 3, 4, 5], "test_data.h5")
        @test true   # If no error occurs, this test passes
    catch e
        @test false  # If an error occurs, this test fails
    end

end
