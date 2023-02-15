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
        directory("resources/STARsolo"),
    threads: 8
    conda:
        "../envs/star.yaml"
    params:
        genomeSAindexNbases=11 if istest else 14,
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeSAindexNbases {params.genomeSAindexNbases}"


# map and count
# https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
rule STARsolo:
    input:
        r1=lambda wc: runs.loc[runs["run_id"] == wc.run, "r1"],
        r2=lambda wc: runs.loc[runs["run_id"] == wc.run, "r2"],
        idx=rules.STARindex.output,
        whitelist=rules.get_whitelist.output,
    output:
        gene_raw=multiext(
            "{outdir}/map_count/{run}/outsGene/raw/",
            "barcodes.tsv",
            "genes.tsv",
            "matrix.mtx",
        ),
        gene_filtered=multiext(
            "{outdir}/map_count/{run}/outsGene/filtered/",
            "barcodes.tsv",
            "genes.tsv",
            "matrix.mtx",
        ),
        genefull_raw=multiext(
            "{outdir}/map_count/{run}/outsGeneFull/raw/",
            "barcodes.tsv",
            "genes.tsv",
            "matrix.mtx",
        ),
        genefull_filtered=multiext(
            "{outdir}/map_count/{run}/outsGeneFull/filtered/",
            "barcodes.tsv",
            "genes.tsv",
            "matrix.mtx",
        ),
        bam="{outdir}/map_count/{run}/Aligned.sortedByCoord.out.bam",  # if outSAMtype="BAM SortedByCoordinate"
        logs=multiext(
            "{outdir}/map_count/{run}/Log", ".out", ".progress.out", ".final.out"
        ),
    threads: 32
    conda:
        "../envs/star.yaml"
    params:
        soloUMIlen=12,  # for 10x v3 chemistry
        soloFeatures="Gene GeneFull",  # Gene=spliced only, GeneFull=spliced and unspliced, and Velocyto=spliced, unspliced, and ambiguous
        soloCellFilter="EmptyDrops_CR",
        soloMultiMappers="EM",
        outSAMtype="BAM SortedByCoordinate",
        soloOutFileNames="outs genes.tsv barcodes.tsv matrix.mtx",
        extra="--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100",
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.r2} {input.r1} \
            --readFilesCommand zcat \
            --soloOutFileNames {params.soloOutFileNames} \
            --soloType CB_UMI_Simple \
            --soloUMIlen {params.soloUMIlen} \
            --soloCBwhitelist {input.whitelist} \
            --soloFeatures {params.soloFeatures} \
            --soloCellFilter {params.soloCellFilter} \
            --soloOutFormatFeaturesGeneField3 "-" \
            --soloMultiMappers {params.soloMultiMappers} \
            --outSAMtype {params.outSAMtype} \
            --outFileNamePrefix "$(dirname {output.bam})/" \
            {params.extra}

        # index output
        sambamba index -t {threads} {output.bam}

        # cleanup tmpdir
        rm -rf $(dirname {output.bam})/_STARtmp
        """


# remove ambient RNA and filter empty droplets
# https://cellbender.readthedocs.io/
rule CellBender:
    input:
        "{outdir}/map_count/{run}/outs{features}/raw/barcodes.tsv",
    output:
        raw="{outdir}/map_count/{run}/outs{features}/cellbender/feature_bc_matrix.h5",
        pdf="{outdir}/map_count/{run}/outs{features}/cellbender/feature_bc_matrix.pdf",
        filtered="{outdir}/map_count/{run}/outs{features}/cellbender/feature_bc_matrix_filtered.h5",
        barcodes="{outdir}/map_count/{run}/outs{features}/cellbender/feature_bc_matrix_barcodes.csv",
    params:
        expected_cells=20000,
        total_droplets_included=100000,
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:latest"
    shell:
        """
        cellbender remove-background \
            --input $(dirname {input}) \
            --output {output.raw} \
            --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets_included}
        """


# quantify transposable element expression
# https://github.com/bodegalab/irescue
# https://www.biorxiv.org/content/10.1101/2022.09.16.508229v2.full
# TODO: not working yet
rule IRescue:
    input:
        bam=rules.STARsolo.output.bam,
        whitelist=rules.CellBender.output.barcodes
        if config["use_CellBender"]
        else "{outdir}/map_count/{run}/outs{features}/filtered/barcodes.tsv",
    output:
        multiext(
            "{outdir}/map_count/{run}/outs{features}/IRescue/",
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
            --genome {params.genome} \
            --threads {threads} \
            --outdir $(dirname {output[0]}) \
            --whitelist {input.whitelist}
        """
