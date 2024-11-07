module adabmDCA
    # using Pkg
    # Pkg.activate("env")

    

    using Base.Threads
    using StatsBase
    using Logging
    using Random
    using Distances
    using LinearAlgebra
    using Distributions
    using Plots
    using Printf
 

    # utils 
    export id, index_interval
    export linear_division, exponential_division_euler, compute_slope
    export oneHotEncoding, oneHot2Categorical, oneHot2Categorical2D
    export remove_autocorrelation, oneHotFreqs, oneHotFij, oneHotFijFast, oneHotFijSymm, oneHotFijSymmFast, oneHotCij, oneHotCijFast
    export read_fasta, read_fasta2, modify_fasta_headers, compute_weights, set_alphabet, read_graph, read_graph_new, initialize_graph, initialize_graph_couplingwise
    export sample_from_profile, gibbs_sampling_edgewise, gibbs_sampling_couplingwise, gibbs_sampling, metropolis_sampling, metropolis_sampling_couplingwise
    export restore_model, restore_model_new, save_model_new, save_chains_new, old2new_model
    export compute_energy, compute_energy_1thread, compute_energy_from_fasta
    export assign_natural_weights
    export wt2dms, distance_from_MSA, oneHotHammingDistance, sample_sid_from_MSA

    # bmDCA
    export fit_bmDCA

    # eaDCA
    export fit_eaDCA

    # edDCA
    export fit_edDCA

    # sample
    export sample_DCA, importance_sample_DCA, TD_integration

    # contacts
    export compute_Frobenius_norm

    include("source/utils.jl")
    include("source/bmDCAsrc.jl")
    include("source/eaDCAsrc.jl")
    include("source/edDCAsrc.jl")
    include("source/resamplesrc.jl")

end