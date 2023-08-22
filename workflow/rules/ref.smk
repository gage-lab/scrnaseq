from snakemake.remote import HTTP, FTP

HTTP = HTTP.RemoteProvider()
FTP = FTP.RemoteProvider()


# download reference genome fasta and gene annotation from 10x Genomics
rule get_refdata:
    input:
        ref10x=HTTP.remote(
            "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
            static=True,
        ),
        rmsk_out=FTP.remote(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz",
            keep_local=True,
            static=True,
        ),
    output:
        fa="resources/genome.fa",
        gtf="resources/genes.gtf",
        rmsk_out="resources/rmsk.out",
    shadow:
        "minimal"
    shell:
        """
        tar -xf {input.ref10x} --wildcards '*genome.fa' '*genes.gtf'
        mv $(basename {input.ref10x} .tar.gz)/genes/genes.gtf {output.gtf}
        mv $(basename {input.ref10x} .tar.gz)/fasta/genome.fa {output.fa}
        gunzip -c {input.rmsk_out} > {output.rmsk_out}
        """


# download barcode whitelist from 10x Genomics
if config["10x_chemistry"] == "3prime_v3":
    whitelist_input = (
        HTTP.remote(
            "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz",
            static=True,
        ),
    )
elif config["10x_chemistry"] == "3prime_v2" or "5prime":
    whitelist_input = (
        HTTP.remote(
            "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt",
            static=True,
        ),
    )


rule get_whitelist:
    input:
        whitelist_input,
    output:
        "resources/10x_whitelist.txt",
    shell:
        """
        if [[ {input} == *.gz ]]; then
            gunzip -c {input} > {output}
        else
            cp {input} {output}
        fi
        """
