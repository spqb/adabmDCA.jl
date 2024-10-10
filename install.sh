#!/bin/bash
echo "Downloading executable scripts..."
wget -O adabmDCA.sh https://raw.githubusercontent.com/spqb/adabmDCA.jl/refs/heads/main/adabmDCA.sh
wget -O execute.jl https://raw.githubusercontent.com/spqb/adabmDCA.jl/refs/heads/main/execute.jl
chmod +x adabmDCA.sh

# Install ArgParse and adabmDCA.jl from the GitHub repo
echo "Installing julia dependencies..."
julia --eval 'using Pkg; Pkg.add("ArgParse"); Pkg.add(PackageSpec(url="https://github.com/spqb/adabmDCA.jl"))'
