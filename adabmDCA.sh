#!/bin/bash

# Check if the first positional argument is provided
if [ -z "$1" ]; then
  echo "Error: No command provided. Use 'train', 'decimate' or 'sample'."
  exit 1
fi

# Assign the first positional argument to a variable
COMMAND=$1
shift # Remove the first positional argument, so "$@" now contains only the optional arguments

# Map the command to the corresponding script
case "$COMMAND" in
  train)
    if [ -z "$1" ]; then
        echo "Error: No sub-command provided for 'train'. Use 'bmDCA', 'eaDCA' or 'edDCA'."
        exit 1
    fi
    
    SUBCOMMAND=$1
    shift
    case "$SUBCOMMAND" in
      -m) 
        SUBCOMMAND=$1
        shift
        case "$SUBCOMMAND" in
          bmDCA)
            exec="bmDCA"
            SUBCOMMAND=$1
            shift
          ;;
          eaDCA)
            exec="eaDCA"
            SUBCOMMAND=$1
            shift
          ;;
          edDCA)
            exec="edDCA"
            SUBCOMMAND=$1
            shift
          ;;
          *)
            echo "1 Error: Invalid sub-command '$SUBCOMMAND' for 'train'. Use -m followed by 'bmDCA', 'eaDCA' or 'edDCA'."
            exit 1
          ;;
        esac
        ;;
      *)
        exec="bmDCA"
        # echo "Error: Invalid sub-command '$SUBCOMMAND' for 'train'. Use -m followed by 'bmDCA', 'eaDCA' or 'edDCA'."
        # exit 1
          ;;
      esac
    ;;
  sample)
    exec="sample"
    ;;
  importance_sample)
    exec="importance_sample"
    ;;
  energies)
    exec="energies"
    ;;
  DMS)
    exec="DMS"
    ;;
  contacts)
    exec="contacts"
    ;;
  *)
    echo "Error: Invalid command '$COMMAND'. Use 'train', 'decimate' or 'sample'."
    exit 1
    ;;
esac

export JULIA_NUM_THREADS=32
# Run the corresponding Julia script with the remaining optional arguments
julia execute.jl -m "$exec" "$SUBCOMMAND""$@" #& disown  # > test.out 
# export JULIA_NUM_THREADS=32
# Run the corresponding Julia script with the remaining optional arguments, loading a module first
# julia -e 'include("src/adabmDCA.jl"); using adamDCA' $train_script "$SUBCOMMAND" "$@" & disown  # > test.out
