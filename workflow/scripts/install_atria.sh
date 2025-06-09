#!/usr/bin/env bash
OUTPUT_PATH=$1
# Create the tools directory if it doesn't exist
mkdir -p $OUTPUT_PATH
TMP_PATH=$(mktemp -d)

export PATH="$HOME/.juliaup/bin:$PATH"  # Add Julia to PATH temporarily
# Install Julia via juliaup and set the default version
if ! command -v julia &> /dev/null; then
    echo "Julia not found, installing Julia"
    curl -fsSL https://install.julialang.org | sh -s -- -y
    juliaup add 1.8
    juliaup default 1.8
else
    echo "Julia is already installed"
fi

# Clone the Atria repository
echo "Cloning Atria repository"
cd $TMP_PATH && \
git clone https://github.com/cihga39871/Atria.git && \
cd Atria && \
git checkout v4.1.1

# Build Atria using the Julia script
if [ ! -f "build_atria.jl" ]; then
    echo "Error: build_atria.jl not found in Atria repository"
    exit 1
fi
julia build_atria.jl $OUTPUT_PATH

cd $OUTPUT_PATH

# Ensure the atria binary exists and is executable
if [ ! -f "atria-4.1.1/bin/atria" ]; then
    echo "Error: Atria binary not found after build"
    exit 1
fi

echo "Atria successfully built and installed at $atria-4.1.1/bin/atria"