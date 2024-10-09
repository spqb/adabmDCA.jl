using Pkg

# Pkg.activate("DCA_env")
include("src/resamplesrc.jl")

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
            default = "DCA_sample"
            help = "(Defaults to DCA_sample). Path to the folder where to save the model."
        "-w", "--weights"
            arg_type = String
            default = nothing
            help = "(Defaults to None). Path to the file containing the weights of the sequences. If None, the weights are computed automatically."
        "-p", "--path_params"
            arg_type = String
            required = true
            help = "(Required) Path to the file containing the model's parameters."
        "-l", "--label"
            arg_type = String
            default = nothing
            help = "label"
        "--alphabet"
            arg_type = String
            default = "protein"
            help = "(Defaults to protein). Type of encoding for the sequences. Choose among ['protein', 'rna', 'dna'] or a user-defined string of tokens."
        "--nsweeps"
            arg_type = Int
            default = 1
            help = "(Defaults to 1). Number of Gibbs sweep for each epoch."
        "--sampler"
            arg_type = String
            default = "gibbs"
            help = "Sampling method to be used. Possible options are gibbs and metropolis"
        "--nchains"
            arg_type = Int
            default = 5_000
            help = "Number of Markov chains to run in parallel. If None, the effective sample size Meff is used with a maximum of 5000"
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
    end

    return parse_args(s)
end

args = parse_commandline()
Threads.nthreads() = args["nthreads"]
println("used threads: ", Threads.nthreads(), "\n")

function main(args)
    println("Parsed args:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end
    println("\n"); flush(stdout)
    resample.resampling_final(args["data"], args["alphabet"], args["weights"], args["nchains"], args["pseudocount"], args["nepochs"], args["nsweeps"], args["output"],  args["path_params"], args["nmeasure"], args["nmix"], args["no_mixingtime"], args["label"], args["showplot"], args["seed"], args["sampler"])
    # resample.resampling_old(args["data"], args["alphabet"], args["weights"], args["nchains"], args["pseudocount"], args["nepochs"], args["nsweeps"], args["output"], args["target"], args["path_params"], args["path_chains"], args["mixingtime"],args["label"], args["sampler"])
end

main(args)