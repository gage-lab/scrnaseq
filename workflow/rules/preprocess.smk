# define filter input based on config values
def get_filter_input(wildcards):
    i = dict()
    if config["use_CellBender"]:
        i["CellBender"] = expand(
            rules.CellBender.output.filtered,
            run=runs["run_id"],
            allow_missing=True,
        )
    else:
        i["STARsolo"] = expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/filtered/matrix.mtx",
            run=runs["run_id"],
            allow_missing=True,
        )
    if config["use_IRescue"]:
        i["IRescue"] = expand(
            "{outdir}/map_count/{run}/IRescue/matrix.mtx.gz",
            run=runs["run_id"],
            allow_missing=True,
        )
    return i


# remove empty droplets, low quality + dead cells, and multiplets
rule filter:
    input:
        unpack(get_filter_input),
    output:
        h5ad="{outdir}/preprocess/filter/{soloFeatures}.h5ad",
        report="{outdir}/preprocess/filter/{soloFeatures}_report.ipynb",
    conda:
        "../envs/pegasus.yaml"
    shadow:
        "shallow"
    log:
        notebook="{outdir}/preprocess/filter/{soloFeatures}_report.ipynb",
    notebook:
        "../notebooks/filter.py.ipynb"


rule render_filter_report:
    input:
        rules.filter.output.report,
    output:
        "{outdir}/preprocess/filter/{soloFeatures}_report.html",
    conda:
        "../envs/jupyter.yaml"
    shell:
        "jupyter nbconvert --no-input --to html {input}"


# rule normalize_pca_scanorama:
#     input:
#         expand(
#             rules.filter.output.mdata,
#             run=runs["run_id"],
#             allow_missing=True,
#         ),
#     output:
#         "{outdir}/integrate/{soloFeatures}.h5ad",
#     conda:
#         "../envs/scanpy.yaml"
#     log:
#         notebook="{outdir}/integrate/{soloFeatures}_normalize_pca_scanorama_report.ipynb",
#     notebook:
#         "../notebooks/normalize_pca_scanorama.py.ipynb"
