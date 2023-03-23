rule install_cplex:
    input:
        directory(config["Compass"]["cplex"]),
    output:
        "resources/cplex_install.done",
    conda:
        "../envs/compass.yaml"
    shell:
        """
        python {input}/python/setup.py install
        touch {output}
        """


# specify output file in directory
rule prepare_data_for_compass:
    input:
        rules.filter.output.h5ad,
    output:
        norm="{outdir}/compass/{soloFeatures}/GRCh38-rna/counts.norm.mtx",
        features="{outdir}/compass/{soloFeatures}/GRCh38-rna/features.tsv",
        barcodes="{outdir}/compass/{soloFeatures}/GRCh38-rna/barcodes.tsv",
        counts="{outdir}/compass/{soloFeatures}/GRCh38-rna/matrix.mtx",
    conda:
        "../envs/pegasus.yaml"
    log:
        "{outdir}/compass/{soloFeatures}/prepare_data_for_compass.log",
    script:
        "../scripts/prepare_data_for_compass.py"


# https://yoseflab.github.io/Compass/tutorial.html
rule run_compass:
    input:
        norm=rules.prepare_data_for_compass.output.norm,
        features=rules.prepare_data_for_compass.output.features,
        barcodes=rules.prepare_data_for_compass.output.barcodes,
        cplex=rules.install_cplex.output,
    output:
        "{outdir}/compass/{soloFeatures}/reactions.tsv.mtx",
    log:
        "{outdir}/compass/{soloFeatures}/compass.log",
    conda:
        "../envs/compass.yaml"
    threads: 8
    shell:
        """
        tmpdir=$(mktemp -d -p $(dirname {output}))
        compass \
            --data-mtx {input.norm} {input.features} {input.barcodes} \
            --output-dir $(dirname {output}) \
            --microcluster-size 10 \
            --num-processes {threads} \
            --temp-dir $tmpdir \
            --species homo_sapiens 2> {log}
        """


rule compass:
    input:
        expand(
            rules.run_compass.output,
            soloFeatures=config["STARsolo"]["soloFeatures"],
            outdir=config["outdir"],
        ),
