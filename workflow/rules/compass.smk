# specify output file in directory
rule prepare_data_for_compass:
    input:
        rules.filter.output.h5ad,
    output:
        directory("{outdir}/compass/{soloFeatures}_cpm.mtx"),
        "{outdir}/compass/{soloFeatures}_cpm.mtx/GRCh38-rna/counts.norm.mtx.gz",
        "{outdir}/compass/{soloFeatures}_cpm.mtx/GRCh38-rna/features.tsv.gz",
        "{outdir}/compass/{soloFeatures}_cpm.mtx/GRCh38-rna/barcodes.tsv.gz",
        "{outdir}/compass/{soloFeatures}_cpm.mtx/GRCh38-rna/matrix.mtx.gz",
    conda:
        "../envs/pegasus.yaml"
    script:
        "../scripts/prepare_data_for_compass.py"


rule run_compass:
    input:
        rules.prepare_data_for_compass.output,
    output:
        "{outdir}/compass/{soloFeatures}/reactions.tsv.mtx",
    conda:
        "../envs/compass.yaml"
    shell:
        "compass --data-mtx {input[1]} {input[2]} {input[3]} --output-dir {wildcards.outdir}/compass/{wildcards.soloFeatures}/ --num-processes 10 --species homo_sapiens"


# https://yoseflab.github.io/Compass/tutorial.html
