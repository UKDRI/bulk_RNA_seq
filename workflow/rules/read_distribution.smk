rule read_distribution:
    input:
        bam="../results/aligned/{sample}.bam",
        annotation=lambda wildcards: f"../results/genome/{selected_genome}/{selected_genome}.bed"
    output:
        txt="../results/qc/rseqc/{sample}_read_distribution.txt"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        set -e  # Exit on first error

        # Ensure output directory exists
        mkdir -p $(dirname {output.txt})

        # Run read_distribution.py
        read_distribution.py -i {input.bam} -r {input.annotation} > {output.txt} 2>&1

        # Verify output exists and is not empty
        if [ ! -s {output.txt} ]; then
            echo "Error: read_distribution.py did not generate expected output." >&2
            exit 1
        fi
        """
