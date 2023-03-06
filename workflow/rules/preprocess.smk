# define filter input based on config values
def get_filter_input(wildcards):
    i = dict()
    if config["use_CellBender"]:
        i["CellBender"] = rules.CellBender.output.filtered
    else:
        i[
            "STARsolo"
        ] = "{outdir}/map_count/{run}/outs{soloFeatures}/filtered/matrix.mtx"
    if config["use_IRescue"]:
        i["IRescue"] = rules.IRescue.output
    return i


# remove empty droplets, lowly expressed genes, low quality + dead cells, and multiplets
rule filter:
    input:
        unpack(get_filter_input),
    output:
        mdata="{outdir}/preprocess/{run}/{soloFeatures}/filtered.h5mu",
        stats="{outdir}/preprocess/{run}/{soloFeatures}/filter_stats.txt",
    conda:
        "../envs/scanpy.yaml"
    log:
        notebook="{outdir}/preprocess/{run}/{soloFeatures}/filter_report.ipynb",
    params:
        expected_multiplet_rate=lambda wc: runs[runs["run_id"] == wc.run][
            "expected_multiplet_rate"
        ].unique()[0],
    notebook:
        "../notebooks/filter.py.ipynb"


rule normalize_pca_scanorama:
    input:
        expand(
            rules.filter.output.mdata,
            run=runs["run_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/integrate/{soloFeatures}.h5ad",
    conda:
        "../envs/scanpy.yaml"
    log:
        notebook="{outdir}/integrate/{soloFeatures}_normalize_pca_scanorama_report.ipynb",
    notebook:
        "../notebooks/normalize_pca_scanorama.py.ipynb"
