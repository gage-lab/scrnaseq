# download reference genome fasta and gene annotation from 10x Genomics
rule get_refdata:
    input:
        ref10x=storage(
            "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
        ),
        rmsk_out=storage(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz",
        ),
    output:
        fa="resources/genome.fa",
        gtf="resources/genes.gtf",
        rmsk_out="resources/rmsk.out",
    shell:
        """
        tar -xf {input.ref10x} --wildcards '*genome.fa' '*genes.gtf'
        mv $(basename {input.ref10x} .tar.gz)/genes/genes.gtf {output.gtf}
        mv $(basename {input.ref10x} .tar.gz)/fasta/genome.fa {output.fa}
        gunzip -c {input.rmsk_out} > {output.rmsk_out}
        """


# download barcode whitelist from Teichlab
if config["10x_chemistry"] == "3prime_v3":
    whitelist_input = (
        "https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz"
    )
elif config["10x_chemistry"] == "3prime_v2" or "5prime":
    whitelist_input = "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz"


rule get_whitelist:
    input:
        storage(whitelist_input),
    output:
        "resources/10x_whitelist.txt",
    shell:
        "gunzip -c {input} > {output}"
