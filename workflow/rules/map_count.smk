# use test reference data if test fastqs are used
if config["istest"]:
    refdata = {
        "fa": "ngs-test-data/scrnaseq_10x_v3/ref/genome.chr21.fa",
        "gtf": "ngs-test-data/scrnaseq_10x_v3/ref/genes.chr21.gtf",
        "rmsk_out": "ngs-test-data/scrnaseq_10x_v3/ref/rmsk_chr21.out",
    }
    assert os.path.exists(refdata["fa"]), "Test reference genome not found"
    assert os.path.exists(refdata["gtf"]), "Test reference annotation not found"
    assert os.path.exists(refdata["rmsk_out"]), "Test repeat masker not found"
else:
    refdata = {
        "fa": rules.get_refdata.output.fa,
        "gtf": rules.get_refdata.output.gtf,
        "rmsk_out": rules.get_refdata.output.rmsk_out,
    }


# reindex genome for faster STAR alignment
rule STARindex:
    input:
        **refdata,
    output:
        protected(sdirectory("resources/STARsolo")),
    threads: 8
    conda:
        "../envs/star.yaml"
    params:
        genomeSAindexNbases=11 if config["istest"] else 14,
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeSAindexNbases {params.genomeSAindexNbases}"


# map and count
# https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# setup STARsolo output
solo_outs = {}
for f in config["STARsolo"]["soloFeatures"]:
    solo_outs[f"{f}_summary"] = "{outdir}/map_count/{run}/outs" + f + "/Summary.csv"
    if f != "Velocyto":
        for d in ["raw", "filtered"]:
            solo_outs[f"{f}_{d}"] = multiext(
                "{outdir}/map_count/{run}/outs" + f + "/" + d + "/",
                "barcodes.tsv",
                "genes.tsv",
                "matrix.mtx",
            )
    else:
        d = "raw"
        solo_outs[f"{f}_{d}"] = multiext(
            "{outdir}/map_count/{run}/outs" + f + "/" + d + "/",
            "barcodes.tsv",
            "genes.tsv",
            "ambiguous.mtx",
            "spliced.mtx",
            "unspliced.mtx",
        )


rule STARsolo:
    input:
        r1=lambda wc: runs.loc[wc.run, "r1"],
        r2=lambda wc: runs.loc[wc.run, "r2"],
        idx=rules.STARindex.output,
        whitelist=rules.get_whitelist.output,
    output:
        **solo_outs,
        bam="{outdir}/map_count/{run}/Aligned.sortedByCoord.out.bam",  # if outSAMtype="BAM SortedByCoordinate"
    threads: 32
    conda:
        "../envs/star.yaml"
    params:
        soloFeatures=" ".join(config["STARsolo"]["soloFeatures"]),
    log:
        multiext("{outdir}/map_count/{run}/Log", ".out", ".progress.out", ".final.out"),
    script:
        "../scripts/STARsolo.py"


rule sambamba_index:
    input:
        rules.STARsolo.output.bam,
    output:
        "{outdir}/map_count/{run}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "{outdir}/map_count/{run}/sambamba_index.log",
    params:
        extra="",  # this must be preset
    threads: 8
    wrapper:
        "v1.23.5/bio/sambamba/index"


rule STARsolo_report:
    input:
        raw=expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/raw/matrix.mtx",
            run=runs["run_id"].unique(),
            allow_missing=True,
        ),
        filtered=expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/filtered/matrix.mtx",
            run=runs["run_id"].unique(),
            allow_missing=True,
        ),
        summary=expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/Summary.csv",
            run=runs["run_id"].unique(),
            allow_missing=True,
        ),
    output:
        "{outdir}/map_count/{soloFeatures}_report.ipynb",
    log:
        notebook="{outdir}/map_count/{soloFeatures}_report.ipynb",
    conda:
        "../envs/pegasus.lock.yaml"
    threads: 1e3
    notebook:
        "../notebooks/STARsolo_report.py.ipynb"


rule render_STARsolo_report:
    input:
        rules.STARsolo_report.output,
    output:
        "{outdir}/map_count/{soloFeatures}_report.html",
    conda:
        "../envs/jupyter.yaml"
    shell:
        "jupyter nbconvert --to html {input}"


# remove ambient RNA and filter empty droplets
# https://cellbender.readthedocs.io/
# check if gpu is available
rule CellBender:
    input:
        raw="{outdir}/map_count/{run}/outs{soloFeatures}/raw/barcodes.tsv",
        filtered="{outdir}/map_count/{run}/outs{soloFeatures}/filtered/barcodes.tsv",
    output:
        raw="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix.h5",
        filtered="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix_filtered.h5",
        pdf="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix.pdf",
        barcodes="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix_cell_barcodes.csv",
    params:
        total_droplets_included=lambda wc: runs[runs["run_id"] == wc.run][
            "total_barcodes"
        ].unique()[0],
    conda:
        # use CUDA if GPU is present
        (
            "../envs/cellbender_cuda.yaml"
            if shutil.which("nvidia-smi")
            else "../envs/cellbender.yaml"
        )
    shell:
        """
        # use CUDA if GPU is present
        if which nvidia-smi > /dev/null; then
            cuda="--cuda"
        else
            cuda=""
        fi
        expected=$(wc -l {input.filtered} | tr ' ' '\n' | head -n 1)

        # run cellbender
        mkdir -p $(dirname {output.raw})
        cellbender remove-background \
            --input $(dirname {input.raw}) \
            --output {output.raw} \
            --epochs 300 \
            --learning-rate 0.00001 \
            --z-dim 50 \
            --expected-cells $expected \
            --total-droplets-included {params.total_droplets_included} $cuda
        """


# quantify transposable element expression
# https://github.com/bodegalab/irescue
# https://www.biorxiv.org/content/10.1101/2022.09.16.508229v2.full
def get_irescue_whitelist(wildcards):
    key = "GeneFull" if "GeneFull" in config["STARsolo"]["soloFeatures"] else "Gene"
    if config["use_CellBender"]:
        return expand(
            rules.CellBender.output.barcodes, soloFeatures=key, allow_missing=True
        )
    else:
        return f"{{outdir}}/map_count/{{run}}/outs{key}/filtered/barcodes.tsv"


rule IRescue:
    input:
        bam=rules.STARsolo.output.bam,
        bai=rules.sambamba_index.output,
        whitelist=get_irescue_whitelist,
    output:
        multiext(
            "{outdir}/map_count/{run}/IRescue/counts/",
            "barcodes.tsv.gz",
            "features.tsv.gz",
            "matrix.mtx.gz",
        ),
    params:
        genome="hg38",
    threads: 8
    log:
        "{outdir}/map_count/{run}/IRescue.log",
    conda:
        "../envs/irescue.yaml"
    shell:
        """
        irescue \
            --bam {input.bam} \
            --genome {params.genome} \
            --threads {threads} \
            --outdir $(dirname $(dirname {output[0]})) \
            --whitelist {input.whitelist} > {log} 2>&1
        """


rule clone_soloTE:
    output:
        directory("resources/soloTE"),
    conda:
        "../envs/solote.yaml"
    log:
        "resources/clone_soloTE.log",
    shell:
        "git clone https://github.com/bvaldebenitom/SoloTE.git {output} 2> {log}"


rule soloTE_build:
    input:
        rules.clone_soloTE.output,
    conda:
        "../envs/solote.yaml"
    output:
        "resources/soloTE_annotation.bed",
    log:
        "resources/soloTE_build.log",
    shadow:
        "shallow"
    shell:
        """
        python {input}/SoloTE_RepeatMasker_to_BED.py -g hg38 && mv hg38_rmsk.bed {output} 2> {log}
        """


rule soloTE:
    input:
        bam=rules.STARsolo.output.bam,
        bai=rules.sambamba_index.output,
        te_bed=rules.soloTE_build.output,
        solote=rules.clone_soloTE.output,
    output:
        multiext(
            "{outdir}/soloTE/{run}_SoloTE_output/",
            "barcodes.tsv",
            "features.tsv",
            "matrix.mtx",
            "{run}_SoloTE.stats",
        ),
    conda:
        "../envs/solote.yaml"
    log:
        "{outdir}/soloTE/{run}_SoloTE_output/soloTE.log",
    threads: 8
    shell:
        """
        outdir=$(dirname $(dirname {output[0]}))

        python {input.solote}/SoloTE_pipeline.py \
            --threads {threads} \
            --bam {input.bam} \
            --teannotation {input.te_bed} \
            --outputprefix {wildcards.run} \
            --outputdir $outdir/ > {log} 2>&1

        mv {wildcards.run}_SoloTE.stats $outdir/{wildcards.run}_SoloTE_output/

        # clean up
        rm -rf $outdir/{wildcards.run}_SoloTE_temp
        """
