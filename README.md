Here's a draft of your README in English for your GitHub package in Julia:

---

# adabmDCA 2.0 - Direct Coupling Analysis in Julia

**Authors:**  
Lorenzo Rosset (Sorbonne Université), Roberto Netti (Sorbonne Université), Anna Paola Muntoni (Politecnico di Torino), Martin Weigt (Sorbonne Université), Francesco Zamponi (Sapienza Università di Roma)  
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

--- 

Let me know if you'd like to make any changes!










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


