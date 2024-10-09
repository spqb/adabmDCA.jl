using Pkg

# Pkg.activate("DCA_env")
include("source/utils.jl")
using .utils
using ArgParse
using Base.Threads

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--path_params"
            arg_type = String
            required = true
            help = "Path to the file containing the model's parameters."
        "-d", "--data"
            arg_type = String
            required = true
            help = "Filename of the dataset to be used for coputing energies using the model."
        "-o", "--output"
            arg_type = String
            default = "sequence_energies.fasta"
            help = "(Defaults to sequence_energies.fasta). Path to the file to save."
        "--nthreads"
            arg_type = Int
            default = 1
            help = "(Defaults to 1). Number of threads used."
    end

    return parse_args(s)
end

# args = parse_commandline()
# Threads.nthreads() = args["nthreads"]
# println("used threads: ", Threads.nthreads(), "\n")

function main(args)
    println("Parsed args:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end
    println("\n"); flush(stdout)
    compute_DMS_energies(args["path_params"], args["data"], args["output"]) 
end

# main(args)