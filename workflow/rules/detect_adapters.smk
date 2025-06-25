rule detect_adapters:
    input:
        fq1=lookup(query="sample_id == '{sample}'", cols="fastq_1", within=samplesheet_df),
        fq2=lookup(query="sample_id == '{sample}'", cols="fastq_2", within=samplesheet_df)
    output:
        adapter_txt="results/adapters/{sample}_adapters.txt"
    log:
        "logs/detect_adapters_{sample}.log"
    container:
        config["atria_image"]
    localrule: 
        not config["trim_reads"]
    params:
        do_trimming = config["trim_reads"]
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.adapter_txt})
        if [[ params.do_trimming == "true" ]]; then

            # Run Atria to detect adapters.
            atria -r {input.fq1} -R {input.fq2} \
                -o $(dirname {output.adapter_txt})/temp_{wildcards.sample} --detect-adapter --threads {threads} \
                > {log} 2>&1

            # Locate the detected adapter summary file
            temp_file=$(find $(dirname {output.adapter_txt})/temp_{wildcards.sample}/ -type f -name "*.txt" | head -n 1)

            # Move the detected adapter file or create an empty one
            if [[ -n "$temp_file" && -f "$temp_file" ]]; then
                mv "$temp_file" {output.adapter_txt}
            else
                echo "Warning: No adapters detected for sample {wildcards.sample}. Creating an empty file." >> {log}
                touch {output.adapter_txt}
            fi
            # Remove the temporary directory
            rm -rf $(dirname {output.adapter_txt})/temp_{wildcards.sample}
        else
            echo "Skipping as 'trim_reads' is false. Creating an empty file." >> {log}
            touch {output.adapter_txt}
        fi
        """