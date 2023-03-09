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


rule integrate:
    input:
        rules.filter.output.h5ad,
    output:
        h5ad="{outdir}/preprocess/integrate/{soloFeatures}.h5ad",
        notebook="{outdir}/preprocess/integrate/{soloFeatures}_report.ipynb",
    threads: 8
    conda:
        "../envs/pegasus.yaml"
    log:
        notebook="{outdir}/preprocess/integrate/{soloFeatures}_report.ipynb",
    notebook:
        "../notebooks/integrate.py.ipynb"


rule render_integrate_report:
    input:
        rules.integrate.output.notebook,
    output:
        "{outdir}/preprocess/integrate/{soloFeatures}_report.html",
    conda:
        "../envs/jupyter.yaml"
    shell:
        "jupyter nbconvert --to html {input}"


rule cluster:
    input:
        rules.integrate.output.h5ad,
    output:
        h5ad="{outdir}/preprocess/cluster/{soloFeatures}.h5ad",
        notebook="{outdir}/preprocess/cluster/{soloFeatures}_report.ipynb",
    threads: 8
    conda:
        "../envs/pegasus.yaml"
    log:
        notebook="{outdir}/preprocess/cluster/{soloFeatures}_report.ipynb",
    notebook:
        "../notebooks/cluster.py.ipynb"


rule render_cluster_report:
    input:
        rules.cluster.output.notebook,
    output:
        "{outdir}/preprocess/cluster/{soloFeatures}_report.html",
    conda:
        "../envs/jupyter.yaml"
    shell:
        "jupyter nbconvert --to html {input}"
