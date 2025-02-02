rule multiqc_trimmed:
    input:
        fastqc_r1=expand("../test-dataset/data/fastqc/trimmed/{sample}_1.atria_fastqc.zip", sample=samples),
        fastqc_r2=expand("../test-dataset/data/fastqc/trimmed/{sample}_2.atria_fastqc.zip", sample=samples)
    output:
        "../test-dataset/data/multiqc/trimmed/multiqc_report.html"
    conda:
        "../envs/multiqc_new.yaml"
    shell:
        """
        # Check if all required FastQC files exist
        missing_files=0
        for file in {input.fastqc_r1} {input.fastqc_r2}; do
            if [ ! -f "$file" ]; then
                echo "Missing file: $file"
                missing_files=1
            fi
        done

        # If any files are missing, create a placeholder report and exit
        if [ $missing_files -eq 1 ]; then
            echo "<html><body><p>MultiQC was skipped due to missing FastQC files.</p></body></html>" > {output}
            exit 0
        fi

        # Run MultiQC if all files are present
        multiqc ../test-dataset/data/fastqc/trimmed/ -o $(dirname {output})
        """
