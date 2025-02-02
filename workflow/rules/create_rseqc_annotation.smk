rule create_rseqc_annotation:
    input:
        gtf="data/genome/mouse/mouse_annotation.gtf"
    output:
        bed12="data/genome/mouse/mouse.bed"
    conda:
        "../envs/rseqc_env.yaml"
    shell:
        """
        # Convert GTF to genePred format
        gtfToGenePred {input.gtf} temp.genePred
        
        # Convert genePred to BED12 format
        genePredToBed temp.genePred {output.bed12}
        
        # Clean up intermediate files
        rm -f temp.genePred
        """
