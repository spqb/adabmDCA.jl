# adabmDCA 2.0 - Direct Coupling Analysis in Julia

**Authors:**  
- **Lorenzo Rosset** (Sorbonne Université, Sapienza Università di Roma)
- **Roberto Netti** (Sorbonne Université)
- **Anna Paola Muntoni** (Politecnico di Torino)
- **Martin Weigt** (Sorbonne Université)
- **Francesco Zamponi** (Sapienza Università di Roma)
  
**Maintainer:** Roberto Netti

## Overview

**adabmDCA 2.0** is a flexible yet easy-to-use implementation of Direct Coupling Analysis (DCA) based on Boltzmann machine learning. This package provides tools for analyzing residue-residue contacts, predicting mutational effects, scoring sequence libraries, and generating artificial sequences, applicable to both protein and RNA families. The package is designed for flexibility and performance, supporting multiple programming languages (C++, Julia, Python) and architectures (single-core/multi-core CPUs and GPUs).  
This repository contains the Julia version of adabmDCA, maintained by **Roberto Netti**.

## Features

- **Direct Coupling Analysis (DCA)** based on Boltzmann machine learning.
- Support for **dense** and **sparse** generative DCA models.
- Available on multiple architectures: single-core and multi-core CPUs, GPUs.
- Ready-to-use for **residue-residue contact prediction**, **mutational-effect prediction**, and **sequence design**.
- Compatible with protein and RNA family analysis.

## Installation

After installing [Julia](https://julialang.org/downloads/) on your system, you can install the package in one of the following ways:

### Option 1: Using the `install.sh` script

1. Download the `install.sh` script to your desired folder.
2. Open a terminal in that folder and run the following commands:

   ```bash
   chmod +x install.sh
   ./install.sh
   ```

This will install all necessary dependencies and set up the package.

### Option 2: Manual Installation via Julia

1. Open Julia and install the package by running:

   ```julia
   using Pkg
   Pkg.add("https://github.com/spqb/adabmDCA.jl")
   Pkg.add("ArgParse")
   ```

2. Download the files `adabmDCA.sh` and `execute.jl` into the same folder.
3. Make the script executable by opening a terminal in the folder and running:

   ```bash
   chmod +x adabmDCA.sh
   ```

This will set up the package for use.

## Usage

To get started with adabmDCA in Julia, please refer to the [documentation](https://github.com/spqb/adabmDCA.jl/wiki) or the [tutorials](https://github.com/spqb/adabmDCA.jl/tree/main/tutorials) for detailed examples of how to apply DCA to real-world protein and RNA datasets.

## License

This package is open-sourced under the MIT License.

## Citation

If you use this package in your research, please cite:

> Rosset, L., Netti, R., Muntoni, A.P., Weigt, M., & Zamponi, F. (2024). adabmDCA 2.0: A flexible but easy-to-use package for Direct Coupling Analysis.

## Acknowledgments

This work was developed in collaboration with Sorbonne Université, Sapienza Università di Roma, and Politecnico di Torino.

