rule process_sj_out:
    input:
        sj_out="../test-dataset/data/aligned_CM/{sample}_SJ.out.tab"
    output:
        splicing_metrics="../test-dataset/data/qc_CM/splicing/{sample}_splicing_metrics.txt"
    conda:
        "../envs/align.yaml"
    shell:
        """
        mkdir -p $(dirname {output.splicing_metrics})

        # Generate summary metrics for splicing analysis
        awk '
            BEGIN {{annotated=0; novel=0; unique=0; multi=0; total=0;}}
            {{
                total++;
                if ($6 == 1) annotated++;
                else novel++;
                unique+=$7;
                multi+=$8;
            }}
            END {{
                print "Sample\\tTotal_Junctions\\tAnnotated_Junctions\\tNovel_Junctions\\tUnique_Reads\\tMulti_Mapped_Reads";
                print "{wildcards.sample}\\t" total "\\t" annotated "\\t" novel "\\t" unique "\\t" multi;
            }}
        ' {input.sj_out} > {output.splicing_metrics}
        """
