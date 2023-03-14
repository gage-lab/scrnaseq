rule install_popscle_helper_tools:
    output:
        directory("resources/popscle_helper_tools"),
    shell:
        """
        git clone https://github.com/aertslab/popscle_helper_tools.git {output}
        """


if "GeneFull" not in config["STARsolo"]["soloFeatures"]:
    demux_barcodes = "{outdir}/map_count/{run}/outsGene/filtered/barcodes.tsv"
    filtered_mtx = ("{outdir}/map_count/{run}/outsGene/filtered/matrix.mtx",)
else:
    demux_barcodes = "{outdir}/map_count/{run}/outsGeneFull/filtered/barcodes.tsv"
    filtered_mtx = ("{outdir}/map_count/{run}/outsGeneFull/filtered/matrix.mtx",)


rule prep_vcf:
    input:
        tools=rules.install_popscle_helper_tools.output,
        bam=expand(
            rules.STARsolo.output.bam, run=runs["run_id"].unique(), allow_missing=True
        ),
        vcf=config["genotypes"],
    output:
        "{outdir}/demuxlet/genotypes.vcf",
    conda:
        "../envs/picard.yaml"
    log:
        "{outdir}/demuxlet/prep_vcf.log",
    shell:
        """
        touch {log} && exec > {log} 2>&1
        source {input.tools}/filter_vcf_file_for_popscle.sh
        check_if_programs_exists

        {input.tools}/sort_vcf_same_as_bam.sh {input.bam[0]} {input.vcf} \
        | only_keep_snps \
        | calculate_AF_AC_AN_values_based_on_genotype_info \
        | filter_out_mutations_missing_genotype_for_one_or_more_samples \
        | filter_out_mutations_homozygous_reference_in_all_samples \
        | filter_out_mutations_homozygous_in_all_samples > {output}
        """


rule filter_for_dsc_pileup:
    input:
        tools=rules.install_popscle_helper_tools.output,
        bam=rules.STARsolo.output.bam,
        barcodes=demux_barcodes,
        vcf=rules.prep_vcf.output,
    output:
        "{outdir}/demuxlet/{run}/filter_for_dsc_pileup.bam",
    conda:
        "../envs/demuxlet.yaml"
    log:
        "{outdir}/demuxlet/{run}/filter_for_dsc_pileup.log",
    threads: 8  # this rule will always 8 threads
    shell:
        "{input.tools}/filter_bam_file_for_popscle_dsc_pileup.sh {input.bam} {input.barcodes} {input.vcf} {output} > {log} 2>&1"


rule dsc_pileup:
    input:
        bam=rules.filter_for_dsc_pileup.output,
        barcodes=demux_barcodes,
        vcf=rules.prep_vcf.output,
    output:
        "{outdir}/demuxlet/{run}/demultiplex.pileup.plp.gz",
    conda:
        "../envs/demuxlet.yaml"
    log:
        "{outdir}/demuxlet/{run}/dsc_pileup.log",
    shell:
        """
        touch {log} && exec > {log} 2>&1

        out=$(dirname {output})/$(basename {output} .plp.gz)

        popscle dsc-pileup \
            --sam {input.bam} \
            --group-list {input.barcodes} \
            --vcf {input.vcf} \
            --out $out
        """


rule demuxlet:
    input:
        plp=rules.dsc_pileup.output,
        barcodes=demux_barcodes,
        vcf=rules.prep_vcf.output,
    output:
        "{outdir}/demuxlet/{run}/demutiplex.demuxlet.best",
    conda:
        "../envs/demuxlet.yaml"
    log:
        "{outdir}/demuxlet/{run}/demuxlet.log",
    shell:
        """
        touch {log} && exec > {log} 2>&1

        out=$(dirname {output})/$(basename {output} .best)
        plp=$(dirname {input.plp})/$(basename {input.plp} .plp.gz)

        popscle demuxlet \
            --plp $plp \
            --group-list {input.barcodes} \
            --vcf {input.vcf} \
            --out $out
        """


rule demuxlet_report:
    input:
        demuxlet=expand(
            rules.demuxlet.output, run=runs["run_id"].unique(), allow_missing=True
        ),
    output:
        "{outdir}/demuxlet/demuxlet_report.ipynb",
    conda:
        "../envs/pegasus.yaml"
    log:
        notebook="{outdir}/demuxlet/demuxlet_report.ipy",
    notebook:
        "../notebooks/demuxlet_report.py.ipynb"


rule render_demuxlet_report:
    input:
        rules.demuxlet_report.output,
    output:
        "{outdir}/demuxlet/demuxlet_report.html",
    conda:
        "../envs/jupyter.yaml"
    shell:
        "juptyer nbconvert --to html {input} --output {output}"
