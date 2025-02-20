rule detect_adapters:
    input:
        fq1="../rawdata/{sample}_1.fq.gz",
        fq2="../rawdata/{sample}_2.fq.gz",
        atria="../results/tools/Atria/atria-4.1.1/bin/atria"
    output:
        adapter_txt="../results/adapters/{sample}_adapters.txt"
    threads: 8
    conda:
        "../envs/qc.yaml"
    shell:
        """
        set -euo pipefail

        # Ensure the output directory exists
        mkdir -p $(dirname {output.adapter_txt})

        # Run Atria to detect adapters
        {input.atria} -r {input.fq1} -R {input.fq2} \
            -o $(dirname {output.adapter_txt})/temp_{wildcards.sample} --detect-adapter --threads {threads}

        # Locate the detected adapter summary file
        temp_file=$(find $(dirname {output.adapter_txt})/temp_{wildcards.sample}/ -type f -name "*.txt" | head -n 1)

        # Move the detected adapter file or create an empty one
        if [[ -n "$temp_file" && -f "$temp_file" ]]; then
            mv "$temp_file" {output.adapter_txt}
            # Remove the temporary directory
            rm -rf $(dirname {output.adapter_txt})/temp_{wildcards.sample}
        else
            echo "Warning: No adapters detected for sample {wildcards.sample}. Creating an empty file."
            touch {output.adapter_txt}
        fi
        """