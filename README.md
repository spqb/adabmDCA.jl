# adabmDCA 2.0 - Direct Coupling Analysis in Julia

## Overview

**adabmDCA 2.0** is a flexible yet easy-to-use implementation of Direct Coupling Analysis (DCA) based on Boltzmann machine learning. This package provides tools for analyzing residue-residue contacts, predicting mutational effects, scoring sequence libraries, and generating artificial sequences, applicable to both protein and RNA families. The package is designed for flexibility and performance, supporting multiple programming languages (C++, Julia, Python) and architectures (single-core/multi-core CPUs and GPUs).  
This repository contains the Julia version of adabmDCA.

The project's main repository can be found at [adabmDCA 2.0](https://github.com/spqb/adabmDCA.git).

> [!NOTE]
>   - 📖 Check out our [Documentation](https://spqb.github.io/adabmDCApy/) website if you want to dive into the package's main features
>   - ❓ Read the reference paper [Rosset et al., 2025](https://doi.org/10.1101/2025.01.31.635874) and its previous version [Muntoni et al., 2021](https://doi.org/10.1186/s12859-021-04441-9) for a detailed description of the proposed methods
>   - 🌐 Explore the [Colab notebook](https://colab.research.google.com/drive/1l5e1W8pk4cB92JAlBElLzpkEk6Hdjk7B?usp=sharing) with a tutorial on some of the package APIs

**Authors:**  
- **Lorenzo Rosset** (Ecole Normale Supérieure ENS, Sorbonne Université)
- **Roberto Netti** (Sorbonne Université)
- **Anna Paola Muntoni** (Politecnico di Torino)
- **Martin Weigt** (Sorbonne Université)
- **Francesco Zamponi** (Sapienza Università di Roma)
  
**Maintainer:** Roberto Netti

## Features

- **Direct Coupling Analysis (DCA)** based on Boltzmann machine learning.
- Support for **dense** and **sparse** generative DCA models.
- Available on multiple architectures: single-core and multi-core CPUs, GPUs.
- Ready-to-use for **residue-residue contact prediction**, **mutational-effect prediction**, and **sequence design**.
- Compatible with protein and RNA family analysis.

## Installation

After installing [Julia](https://julialang.org/downloads/) on your system, you can install the package in one of the following ways:

### Option 1: Using bash command
Open a terminal in the desired folder, and run the following commands:

   ```bash
   # Download scripts from Github
   wget -O adabmDCA.sh https://raw.githubusercontent.com/spqb/adabmDCA.jl/refs/heads/main/adabmDCA.sh
   wget -O execute.jl https://raw.githubusercontent.com/spqb/adabmDCA.jl/refs/heads/main/execute.jl
   chmod +x adabmDCA.sh

   # Install ArgParse and adabmDCA.jl from the GitHub repo
   julia --eval 'using Pkg; Pkg.add("ArgParse"); Pkg.add(PackageSpec(url="https://github.com/spqb/adabmDCA.jl"))'
   ```
This will install all necessary dependencies and set up the package.

Here’s the refined text for your README:

---

If you want to use the simpler `adabmDCA` command instead of typing `./adabmDCA.sh`, you can create an alias in your shell configuration file. Add the following line to your `~/.bashrc` (or `~/.zshrc` for Zsh):

```bash
alias adabmDCA="$(pwd)/adabmDCA.sh"
```

Then, reload your shell configuration:

```bash
source ~/.bashrc  # Use this for Bash
# source ~/.zshrc  # Use this for Zsh
```

**Note:** You need to be in the same directory as the `adabmDCA.sh` and `execute.jl` files when creating the alias, as the command uses the current folder path.

---

Feel free to copy and paste it into your README! 😊


### Option 2: Manual Installation via Julia

1. Open Julia and install the package by running:

   ```julia
   using Pkg
   Pkg.add(url="https://github.com/spqb/adabmDCA.jl")
   Pkg.add("ArgParse")
   ```

2. Download the files `adabmDCA.sh` and `execute.jl` into the same folder.
3. Make the script executable by opening a terminal in the folder and running:

   ```bash
   chmod +x adabmDCA.sh
   ```

This will set up the package for use.

## Usage

To get started with adabmDCA in Julia, please refer to the [documentation](https://spqb.github.io/adabmDCApy/).
## License

This package is open-sourced under the MIT License.

## Citation

If you use this package in your research, please cite:

> Rosset, L., Netti, R., Muntoni, A.P., Weigt, M., & Zamponi, F. (2024). adabmDCA 2.0: A flexible but easy-to-use package for Direct Coupling Analysis.


## Acknowledgments

This work was developed in collaboration with Sorbonne Université, Sapienza Università di Roma, and Politecnico di Torino.

