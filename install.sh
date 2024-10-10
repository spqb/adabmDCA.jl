#!/bin/bash

wget -O adabmDCA.sh

chmod +x adabmDCA.sh

# Install ArgParse and adabmDCA.jl from the GitHub repo
julia --eval 'using Pkg; Pkg.add("ArgParse"); Pkg.add(PackageSpec(url="https://github.com/spqb/adabmDCA.jl"))'
