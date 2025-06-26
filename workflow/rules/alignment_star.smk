adjust_ulimit()

rule alignment_star:
    """
    Align paired-end reads with STAR for a sample that includes {sample} and {genome}.
    """
    input:
        # Trimmed FASTQs
        r1 = "results/trimmed/{sample}_1.atria.fastq.gz",
        r2 = "results/trimmed/{sample}_2.atria.fastq.gz",
        # Use a lambda to derive the star index path from the genome wildcard
        star_index=expand(config["star_index"], selected_genome=config["selected_genome"])
    output:
        bam = ensure("results/aligned/{sample}.bam", non_empty=True),
        log = ensure("results/aligned/{sample}_Log.out", non_empty=True),
        sj  = ensure("results/aligned/{sample}_SJ.out.tab", non_empty=True)
    log:
        "logs/alignment_star_{sample}.log"
    conda:
        "../envs/align.yaml"
    params:
        tmp_dir = get_star_tempdir,
        star_ram = config["star_ram_limit"] # Limit for STAR's RAM usage in bytes
    threads: 10
    resources:
        mem_gb = config["star_ram_limit"] / 1000000000
        runtime_h = 1
    shell:
        """
        STAR \
          --genomeDir "{input.star_index}" \
          --readFilesIn "{input.r1}" "{input.r2}" \
          --readFilesCommand zcat \
          --runThreadN {threads} \
          --outFileNamePrefix "{params.tmp_dir}/" \
          --outSAMtype BAM SortedByCoordinate \
          --limitBAMsortRAM {params.star_ram} \
          > {log} 2>&1

        mv "{params.tmp_dir}/Aligned.sortedByCoord.out.bam" "{output.bam}"
        mv "{params.tmp_dir}/Log.final.out" "{output.log}"
        mv "{params.tmp_dir}/SJ.out.tab" "{output.sj}"

        rm -rf "{params.tmp_dir}"
        """
