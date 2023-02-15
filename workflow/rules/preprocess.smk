# define filter input based on config values
def get_filter_input(wildcards):
    i = dict()
    if config["use_CellBender"]:
        i["CellBender"] = rules.CellBender.output.filtered
    else:
        i["STARsolo"] = "{outdir}/map_count/{run}/outs{features}/filtered/matrix.mtx"
    if config["use_IRescue"]:
        i["IRescue"] = rules.IRescue.output
    return i


# remove empty droplets, lowly expressed genes, low quality + dead cells, and multiplets
rule filter:
    input:
        unpack(get_filter_input),
    output:
        "{outdir}/preprocess/{run}/{features}/filtered.h5mu",
    conda:
        "../envs/scanpy.yaml"
    log:
        notebook="{outdir}/preprocess/{run}/{features}/filter_report.ipynb",
    notebook:
        "../notebooks/filter.py.ipynb"
