rule download_atria:
    output:
        "../results/tools/Atria/atria-4.1.1/bin/atria"
    conda:
        "../envs/trim-env.yaml"  # Correct environment path
    shell:
        """
        # Create the tools directory if it doesn't exist
        mkdir -p tools

        cd tools

        # Install Julia via juliaup and set the default version
        if ! command -v julia &> /dev/null; then
            echo "Julia not found, installing Julia"
            curl -fsSL https://install.julialang.org | sh
            export PATH="$HOME/.juliaup/bin:$PATH"  # Add Julia to PATH temporarily
            juliaup add 1.8
            juliaup default 1.8
        else
            echo "Julia is already installed"
        fi

        # Remove existing Atria directory if it exists
        if [ -d "Atria" ]; then
            echo "Atria directory exists, removing it"
            rm -rf Atria
        fi

        # Clone the Atria repository
        echo "Cloning Atria repository"
        git clone https://github.com/cihga39871/Atria.git
        cd Atria

        # Build Atria using the Julia script
        if [ ! -f "build_atria.jl" ]; then
            echo "Error: build_atria.jl not found in Atria repository"
            exit 1
        fi
        julia build_atria.jl

        # Ensure the atria binary exists and is executable
        if [ ! -f "atria-4.1.1/bin/atria" ]; then
            echo "Error: Atria binary not found after build"
            exit 1
        fi
        chmod +x atria-4.1.1/bin/atria

        echo "Atria successfully built and installed at ../results/tools/Atria/atria-4.1.1/bin/atria"
        """
