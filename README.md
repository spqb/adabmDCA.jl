# adabmDCA 2.0

## Overview

**adabmDCA 2.0** is a flexible yet easy-to-use implementation of Direct Coupling Analysis (DCA), designed for analyzing protein and RNA families. The package leverages Boltzmann machine learning and offers support for both dense and sparse generative models. It also facilitates various downstream tasks, including:

- Residue-residue contact prediction
- Mutational-effect prediction
- Scoring of sequence libraries
- Generation of artificial sequences for sequence design

This package is available in multiple programming languages (C++, Julia, Python) and can be run on different architectures (single-core/multi-core CPUs, GPUs), providing a common front-end interface for users. 

Developed by:
- **Lorenzo Rosset** (Sorbonne Université, CNRS)
- **Roberto Netti** (Sorbonne Université, CNRS)
- **Anna Paola Muntoni** (Politecnico di Torino)
- **Martin Weigt** (Sorbonne Université, CNRS)
- **Francesco Zamponi** (Sapienza Università di Roma)

## Features

- Flexible implementation of DCA using Boltzmann machine learning
- Supports dense and sparse models
- Usable for protein and RNA families
- Compatible with multiple architectures (CPU/GPU)
- Provides tools for common bioinformatics tasks

## Installation

To get started with **adabmDCA 2.0**, follow these steps:

### Step 1: Install Julia

If you haven't installed Julia yet, follow the official instructions [here](https://julialang.org/downloads/).

### Step 2: Install the Required Dependencies

1. Download the `install.sh` file to your desired directory.
2. Open a terminal in that directory and run the following commands:

   ```bash
   chmod +x install.sh
   ./install.sh
   ```

This script will install all the necessary dependencies.

### Step 3: Install the Package

You can manually install **adabmDCA 2.0** in Julia using the following command:

```julia
using Pkg
Pkg.add("https://github.com/spqb/adabmDCA.jl")
```

After installation, download the `adabmDCA.sh` and `execute.jl` files into the same directory.

### Step 4: Final Setup

1. Open a terminal in the directory where you downloaded `adabmDCA.sh`.
2. Make the script executable:

   ```bash
   chmod +x adabmDCA.sh
   ```

You are now ready to use **adabmDCA 2.0**!

## Usage

A tutorial on how to use **adabmDCA 2.0** is provided within the package. You can explore various learning protocols for different models and apply them to protein or RNA sequence data.

## License

This project is licensed under the MIT License.

## Authors

- Roberto Netti – Sorbonne Université
- Lorenzo Rosset – Sorbonne Université
- Anna Paola Muntoni – Politecnico di Torino
- Martin Weigt – Sorbonne Université
- Francesco Zamponi – Sapienza Università di Roma

## Citation

If you use this package in your research, please cite:

> Rosset, L., Netti, R., Muntoni, A.P., Weigt, M., & Zamponi, F. (2024). adabmDCA 2.0: A flexible but easy-to-use package for Direct Coupling Analysis.

