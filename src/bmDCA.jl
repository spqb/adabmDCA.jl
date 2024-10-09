using Pkg

# Pkg.activate("DCA_env")
include("source/bmDCAsrc.jl")

using .bmDCAsrc
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
        "--graph"
            arg_type = String
            default = nothing
            help = "External graph on which Boltzmann learning is performed."
    end

    return parse_args(s)
end

# args = parse_commandline()
# Threads.nthreads() = args["nthreads"]
# println("used threads: ", Threads.nthreads())

function main(args)
    println("Parsed args:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end
    println("\n"); flush(stdout)
    fit_bmDCA(args["data"], args["alphabet"], args["weights"], args["nchains"], args["pseudocount"], args["lr"], args["nepochs"], args["nsweeps"], args["output"], args["target"], args["graph"], args["path_params"], args["path_chains"], args["label"], args["sampler"], args["seed"])
end

# main(args)