# define filter input based on config values
def get_filter_input(wildcards):
    i = dict()
    if config["use_CellBender"]:
        i["CellBender"] = expand(
            rules.CellBender.output.filtered, run=runs["run_id"], allow_missing=True
        )
    else:
        i["STARsolo"] = expand(
            "{outdir}/map_count/{run}/outs{features}/filtered/matrix.mtx",
            run=runs["run_id"],
            allow_missing=True,
        )
    if config["use_IRescue"]:
        i["IRescue"] = expand(
            rules.IRescue.output, run=runs["run_id"], allow_missing=True
        )
    return i


# remove empty droplets, lowly expressed genes, low quality + dead cells, and multiplets
rule filter:
    input:
        unpack(get_filter_input),
    output:
        "{outdir}/preprocess/{features}/filtered.h5mu",
    conda:
        "../envs/scanpy.yaml"
    log:
        notebook="{outdir}/preprocess/{features}/filter_report.ipynb",
    notebook:
        "../notebooks/filter.py.ipynb"
