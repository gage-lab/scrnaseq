# download barcode whitelist from Teichlab
if config["10x_chemistry"] == "3prime_v3":
    whitelist_input = "https://github.com/Teichlab/scg_lib_structs/raw/refs/heads/master/data/10X-Genomics/3M-february-2018.txt.gz"
elif config["10x_chemistry"] == "3prime_v2" or "5prime":
    whitelist_input = "https://github.com/Teichlab/scg_lib_structs/raw/refs/heads/master/data/10X-Genomics/737K-august-2016.txt.gz"
else:
    ValueError("10x_chemistry must be 3prime_v2, 3prime_v3, or 5prime")


rule get_whitelist:
    input:
        storage(whitelist_input),
    output:
        "resources/10x_whitelist.txt",
    shell:
        "gunzip -c {input} > {output}"
