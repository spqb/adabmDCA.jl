# using Pkg

# Pkg.activate("env")
# include("src/adabmDCA.jl")
using adabmDCA

using ArgParse
# using Base.Threads


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-m", "--model"
            arg_type = String
            required = true
            help = "method."
        "-d", "--data"
            arg_type = String
            required = false
            help = "Filename of the dataset to be used for training the model."
        "-o", "--output"
            arg_type = String
            default = "DCA_model"
            help = "(Defaults to DCA_model). Path to the folder where to save the model."
        "-w", "--weights"
            arg_type = String
            default = nothing
            help = "(Defaults to None). Path to the file containing the weights of the sequences. If None, the weights are computed automatically."
        "-p", "--path_params"
            arg_type = String
            default = nothing
            help = "(Defaults to None) Path to the file containing the model's parameters. Required for restoring the training."
        "-c", "--path_chains"
            arg_type = String
            default = nothing
            help = "(Defaults to None) Path to the fasta file containing the model's chains. Required for restoring the training."
        "-l", "--label"
            arg_type = String
            default = nothing
            help = "label"
        "--alphabet"
            arg_type = String
            default = "protein"
            help = "(Defaults to protein). Type of encoding for the sequences. Choose among ['protein', 'rna', 'dna'] or a user-defined string of tokens."
        "--lr"                 
            arg_type = Float64
            default = 0.01
            help = "(Defaults to 0.01). Learning rate."
        "--nsweeps"
            arg_type = Int
            default = 10
            help = "(Defaults to 10). Number of Gibbs steps for each gradient estimation."
        "--sampler"
            arg_type = String
            default = "gibbs"
            help = "Sampling method to be used. Possible options are gibbs and metropolis"
        "--nchains"
            arg_type = Int
            default = 5_000
            help = "Number of Markov chains to run in parallel. If None, the effective sample size Meff is used with a maximum of 5000"
        "--target"
            arg_type = Float64
            default = 0.95
            help = "(Defaults to 0.95). Pearson correlation coefficient on the two-points statistics to be reached."
        "--nepochs"
            arg_type = Int
            default = 50_000
            help = "(Defaults to 50_000). Maximum number of epochs allowed."
        "--pseudocount"
            arg_type = Float32
            default = nothing
            help = "(Defaults to None). Pseudo count for the single and two points statistics. If None, it is set to 1/Meff."   
        "--nthreads"
            arg_type = Int
            default = 1
            help = "(Defaults to 1). Number of threads used."
        "--seed"
            arg_type = Int
            default = 0
            help = "(Defaults to 0). Random seed."
        "--graph"
            arg_type = String
            default = nothing
            help = "External graph on which Boltzmann learning is performed."

        # eaDCA
        "--gsteps"
            arg_type = Int
            default = 10
            help = "(Defaults to 10). Number of gradient updates between two subsequent coupling activations."
        "--factivate"
            arg_type = Float64
            default = 0.001
            help = "(Defaults to 100). Fraction of coulings activated at each epoch."

        # edDCA
        "--density"
            arg_type = Float64
            default = 0.02
            help = "(Defaults to 0.02). Target density to be reached."
        "--drate"
            arg_type = Float64
            default = 0.01
            help = "(Defaults to 0.01). Fraction of remaining couplings to be pruned at each decimation step."
        "--max_convergence_step"
            arg_type = Int
            default = 10000
            help = "(Defaults to 10000). Maximum number of convergence step."
        "--target_densities"
            arg_type = String
            default = nothing
            help = "(Defaults to nothing). Array containing target densities"

        # sample
        "--no_mixingtime"
            default = true
            help = "(Defaults to True) Compute mixing time."
            action = :store_false
        "--nmeasure"
            arg_type = Int
            default = 5_000
            help = "Number of Markov chains to run in parallel for mixing time estimation."
        "--nmix"
            arg_type = Int
            default = 2
            help = "Sampling will be done for a total time of nmix * t_mix"
        "--showplot"
            default = false
            help = "(Defaults to False) show plot"
            action = :store_true

        # importance sampling
        "--theta"
            arg_type = Float32
            default = 1
            help = "Bias towards target sequence"
        "--targetseq"
            arg_type = String
            help = "Filename of the dataset to be used for training the model."

        # TD integration 
        "--intstep"
            arg_type = Int64
            default = 100
            help = "integration step"
    end

    return parse_args(s)
end

args = parse_commandline()
Threads.nthreads() = args["nthreads"]
println("used threads: ", Threads.nthreads())

if args["model"] == "bmDCA"
    fit_bmDCA(args["data"], args["alphabet"], args["weights"], args["nchains"], args["pseudocount"], args["lr"], args["nepochs"], args["nsweeps"], args["output"], args["target"], args["graph"], args["path_params"], args["path_chains"], args["label"], args["sampler"], args["seed"])

elseif args["model"] == "eaDCA"
    fit_eaDCA(args["data"], args["alphabet"], args["weights"], args["nchains"], args["factivate"],args["pseudocount"], args["lr"], args["nepochs"], args["nsweeps"], args["gsteps"], args["output"], args["target"], args["path_params"], args["path_chains"], args["label"], args["sampler"], args["seed"])

elseif args["model"] == "edDCA"
    fit_edDCA(args["data"], args["path_params"], args["path_chains"], args["density"],  args["target"], args["drate"], args["lr"], args["nsweeps"], args["nchains"], args["alphabet"], args["weights"], args["pseudocount"], args["nepochs"], args["output"], args["max_convergence_step"], args["label"], args["sampler"], args["seed"], args["target_densities"])

elseif args["model"] == "sample"
    sample_DCA(args["data"], args["alphabet"], args["weights"], args["nchains"], args["pseudocount"], args["nepochs"], args["nsweeps"], args["output"],  args["path_params"], args["nmeasure"], args["nmix"], args["no_mixingtime"], args["label"], args["showplot"], args["seed"], args["sampler"])

elseif args["model"] == "energies"
    compute_energy_from_fasta(args["path_params"], args["data"], args["output"], args["label"])

elseif args["model"] == "DMS"
    compute_DMS_energies(args["path_params"], args["data"], args["output"]) 

elseif args["model"] == "contacts"
    compute_Frobenius_norm(args["path_params"], args["output"], args["label"]) 




elseif args["model"] == "importance_sample"
    importance_sample_DCA(args["data"], args["alphabet"], args["weights"], args["nchains"], args["pseudocount"], args["nepochs"], args["nsweeps"], args["output"],  args["path_params"], args["nmeasure"], args["nmix"], args["no_mixingtime"], args["label"], args["showplot"], args["seed"], args["sampler"], args["targetseq"], args["theta"])
elseif args["model"] == "TD_integration"
    TD_integration(args["data"], args["alphabet"], args["weights"], args["nchains"], args["nsweeps"], args["output"],  args["path_params"], args["path_chains"], args["label"], args["sampler"], args["targetseq"], args["intstep"], args["theta"])
end