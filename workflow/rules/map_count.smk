# reindex genome for faster STAR alignment
rule STARindex:
    input:
        **refdata,
    output:
        directory("resources/STARsolo"),
    threads: 8
    conda:
        "../envs/star.yaml"
    params:
        genomeSAindexNbases=11 if config["istest"] else 14,
        gff=(
            "--sjdbGTFtagExonParentTranscript Parent"
            if ".gff" in refdata["gtf"]
            else ""
        ),
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} {params.gff} \
            --genomeSAindexNbases {params.genomeSAindexNbases}
        """


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
