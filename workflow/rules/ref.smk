from snakemake.remote import HTTP

HTTP = HTTP.RemoteProvider()


# download reference genome fasta and gene annotation from 10x Genomics
rule get_refdata:
    input:
        HTTP.remote(
            "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
            static=True,
        ),
    output:
        fa="resources/genome.fa",
        gtf="resources/genes.gtf",
    shadow:
        "minimal"
    shell:
        """
        tar -xf {input} --wildcards '*genome.fa' '*genes.gtf'
        mv $(basename {input} .tar.gz)/genes/genes.gtf {output.gtf}
        mv $(basename {input} .tar.gz)/fasta/genome.fa {output.fa}
        """


# download barcode whitelist from 10x Genomics
# for 10x v3 chemistry
rule get_whitelist:
    input:
        HTTP.remote(
            "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz",
            static=True,
        ),
    output:
        "resources/10x_v3_whitelist.txt",
    shell:
        "gunzip -c {input} > {output}"
