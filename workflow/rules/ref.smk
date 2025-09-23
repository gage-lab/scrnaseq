# download barcode whitelist from Teichlab
if CHEMISTRY == "3prime_v3":
    whitelist_input = "https://github.com/Teichlab/scg_lib_structs/raw/refs/heads/master/data/10X-Genomics/3M-february-2018.txt.gz"
elif CHEMISTRY == "3prime_v2" or CHEMISTRY == "5prime":
    whitelist_input = "https://github.com/Teichlab/scg_lib_structs/raw/refs/heads/master/data/10X-Genomics/737K-august-2016.txt.gz"
elif CHEMISTRY == "3prime_v4":
    whitelist_input = "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-3pgex-may-2023.txt.gz"
else:
    ValueError("10x_chemistry must be 3prime_v2, 3prime_v3, 3prime_v4, or 5prime")


rule get_whitelist:
    input:
        storage(whitelist_input),
    output:
        f"resources/10x_whitelist_{CHEMISTRY}.txt",
    shell:
        "gunzip -c {input} > {output}"
