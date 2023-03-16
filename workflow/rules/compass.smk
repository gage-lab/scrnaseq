rule prepare_data_for_compass:
    input:
        rules.filter.output.h5ad,
    output:
        directory("{outdir}/compass/{soloFeatures}_cpm.mtx"),
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
        "compass --data-mtx {input[0]}/counts.norm.mtx.gz {input[0]}/barcodes.tsv.gz {input[0]}/features.tsv.gz --output-dir {wildcards.outdir}/compass/{wildcards.soloFeatures}/ --num-processes 10 --species homo_sapiens"


# https://yoseflab.github.io/Compass/tutorial.html
