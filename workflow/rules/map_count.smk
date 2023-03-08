# use test reference data if test fastqs are used
if istest:
    refdata = {
        "fa": "ngs-test-data/scrnaseq_10x_v3/ref/genome.chr21.fa",
        "gtf": "ngs-test-data/scrnaseq_10x_v3/ref/genes.chr21.gtf",
    }
    assert os.path.exists(refdata["fa"]), "Test reference genome not found"
    assert os.path.exists(refdata["gtf"]), "Test reference annotation not found"
else:
    refdata = {"fa": rules.get_refdata.output.fa, "gtf": rules.get_refdata.output.gtf}


# reindex genome for faster STAR alignment
rule STARindex:
    input:
        **refdata,
    output:
        protected(directory("resources/STARsolo")),
    threads: 8
    conda:
        "../envs/star.yaml"
    params:
        genomeSAindexNbases=11 if istest else 14,
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeSAindexNbases {params.genomeSAindexNbases}"


# map and count
# https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# setup STARsolo output
solo_outs = {}
for f in config["STARsolo"]["soloFeatures"]:
    solo_outs[f"{f}_summary"] = f"{{outdir}}/map_count/{{run}}/outs{f}/Summary.csv"
    for d in ["raw", "filtered"]:
        if f != "Velocyto":
            solo_outs[f"{f}_{d}"] = multiext(
                f"{{outdir}}/map_count/{{run}}/outs{f}/{d}/",
                "barcodes.tsv",
                "genes.tsv",
                "matrix.mtx",
            )
        else:
            solo_outs[f"{f}_{d}"] = multiext(
                f"{{outdir}}/map_count/{{run}}/outs{f}/{d}/",
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
        logs=multiext(
            "{outdir}/map_count/{run}/Log", ".out", ".progress.out", ".final.out"
        ),
    threads: 32
    conda:
        "../envs/star.yaml"
    params:
        soloFeatures=" ".join(config["STARsolo"]["soloFeatures"]),
    script:
        "../scripts/STARsolo.py"


rule STARsolo_report:
    input:
        raw=expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/raw/matrix.mtx",
            run=runs["run_id"],
            allow_missing=True,
        ),
        filtered=expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/filtered/matrix.mtx",
            run=runs["run_id"],
            allow_missing=True,
        ),
        summary=expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/Summary.csv",
            run=runs["run_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/map_count/{soloFeatures}_report.ipynb",
    log:
        notebook="{outdir}/map_count/{soloFeatures}_report.ipynb",
    conda:
        "../envs/pegasus.yaml"
    notebook:
        "../notebooks/STARsolo_report.py.ipynb"


# remove ambient RNA and filter empty droplets
# https://cellbender.readthedocs.io/
# check if gpu is available
rule CellBender:
    input:
        "{outdir}/map_count/{run}/outs{soloFeatures}/raw/barcodes.tsv",
    output:
        raw="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix.h5",
        filtered="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix_filtered.h5",
        pdf="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix.pdf",
        barcodes="{outdir}/map_count/{run}/outs{soloFeatures}/cellbender/feature_bc_matrix_cell_barcodes.csv",
    params:
        expected_cells=lambda wc: runs[runs["run_id"] == wc.run][
            "expected_cells"
        ].unique()[0],
        total_droplets_included=lambda wc: runs[runs["run_id"] == wc.run][
            "total_barcodes"
        ].unique()[0],
    conda:
        # use CUDA if GPU is present
        "../envs/cellbender_cuda.yaml" if shutil.which(
        "nvidia-smi"
        ) else "../envs/cellbender.yaml"
    shell:
        """
        # use CUDA if GPU is present
        if which nvidia-smi > /dev/null; then
            cuda="--cuda"
        else
            cuda=""
        fi

        # run cellbender
        mkdir -p $(dirname {output.raw})
        cellbender remove-background \
            --input $(dirname {input}) \
            --output {output.raw} \
            --expected-cells {params.expected_cells} \
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
        whitelist=get_irescue_whitelist,
    output:
        multiext(
            "{outdir}/map_count/{run}/IRescue/",
            "barcodes.tsv.gz",
            "features.tsv.gz",
            "matrix.mtx.gz",
        ),
    params:
        genome="hg38",
    threads: 8
    conda:
        "../envs/irescue.yaml"
    shell:
        """
        irescue \
            --bam {input.bam} \
            --tmpdir IRescue_tmp_{wildcards.run} \
            --genome {params.genome} \
            --threads {threads} \
            --outdir $(dirname {output[0]}) \
            --whitelist {input.whitelist}
        """
