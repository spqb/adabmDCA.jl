using Pkg

# Pkg.activate("DCA_env")
include("source/edDCAsrc.jl")

using .edDCAsrc
using ArgParse
using Base.Threads

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "-d", "--data"
            arg_type = String
            required = true
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
    end
    return parse_args(s)
end

args = parse_commandline()
Threads.nthreads() = args["nthreads"]
println("used threads: ", Threads.nthreads(), "\n"); flush(stdout)

function main(args)
    println("Parsed args:"); flush(stdout)
    for (arg,val) in args
        println("  $arg  =>  $val"); flush(stdout)
    end
    println("\n"); flush(stdout)
    fit_edDCA(args["data"], args["path_params"], args["path_chains"], args["density"],  args["target"], args["drate"], args["lr"], args["nsweeps"], args["nchains"], args["alphabet"], args["weights"], args["pseudocount"], args["nepochs"], args["output"], args["max_convergence_step"], args["label"], args["sampler"], args["seed"])
end

main(args)