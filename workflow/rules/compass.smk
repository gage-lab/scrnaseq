rule prepare_data_for_compass:
    input:
        rules.filter.output.h5ad,
    output:
        "{outdir}/compass/{soloFeatures}_cpm.mtx",
    conda:
        "../envs/pegasus.yaml"
    script:
        "scripts/prepare_data_for_compass.py"


# https://yoseflab.github.io/Compass/tutorial.html
