# define filter input based on config values
def get_filter_input(wildcards):
    i = dict()

    # choose the input matrix
    if config["use_CellBender"]:
        i["CellBender"] = expand(
            rules.CellBender.output.filtered,
            run=runs["run_id"].unique(),
            allow_missing=True,
        )
    else:
        i["STARsolo"] = expand(
            "{outdir}/map_count/{run}/outs{soloFeatures}/filtered/matrix.mtx",
            run=runs["run_id"].unique(),
            allow_missing=True,
        )

    # add demuxlet if available
    if config["use_Demuxlet"]:
        # validate demuxlet config
        assert (
            "genotypes" in config.keys()
        ), "Demuxlet is activated by genotype file is not provided"
        assert os.path.isfile(config["genotypes"]), (
            config["genotypes"] + " file does not exist"
        )
        assert config["genotypes"].endswith(".vcf"), (
            config["genotypes"] + " is not a vcf file"
        )

        i["Demuxlet"] = expand(
            rules.demuxlet.output,
            run=runs["run_id"].unique(),
            allow_missing=True,
        )
        i["Demuxlet_report"] = rules.demuxlet_report.output

    return i


# remove empty droplets, low quality + dead cells, and multiplets
rule filter:
    input:
        unpack(get_filter_input),
    output:
        h5ad="{outdir}/preprocess/01_filter/{soloFeatures}.h5ad",
        zarr="{outdir}/preprocess/01_filter/{soloFeatures}.zarr.zip",
        seurat="{outdir}/preprocess/01_filter/{soloFeatures}.rds",
    conda:
        "../envs/pegasus.yaml"
    log:
        notebook="{outdir}/preprocess/01_filter/{soloFeatures}_report.ipynb",
    threads: 1e3  # force only one to run at a time, to prevent tmpfile conflicts
    notebook:
        "../notebooks/filter.py.ipynb"


rule integrate:
    input:
        rules.filter.output.zarr,
    output:
        h5ad="{outdir}/preprocess/02_integrate/{soloFeatures}.h5ad",
        zarr="{outdir}/preprocess/02_integrate/{soloFeatures}.zarr.zip",
        seurat="{outdir}/preprocess/02_integrate/{soloFeatures}.rds",
    threads: 8
    conda:
        "../envs/pegasus.yaml"
    log:
        notebook="{outdir}/preprocess/02_integrate/{soloFeatures}_report.ipynb",
    notebook:
        "../notebooks/integrate.py.ipynb"


rule cluster:
    input:
        rules.integrate.output.zarr,
    output:
        h5ad="{outdir}/preprocess/03_cluster/{soloFeatures}.h5ad",
        zarr="{outdir}/preprocess/03_cluster/{soloFeatures}.zarr.zip",
        seurat="{outdir}/preprocess/03_cluster/{soloFeatures}.rds",
    threads: 8
    conda:
        "../envs/pegasus.yaml"
    log:
        notebook="{outdir}/preprocess/03_cluster/{soloFeatures}_report.ipynb",
    notebook:
        "../notebooks/cluster.py.ipynb"


rule render_preprocess_reports:
    input:
        "{outdir}/preprocess/{stage}/{soloFeatures}_report.ipynb",
    output:
        "{outdir}/preprocess/{stage}/{soloFeatures}_report.html",
    conda:
        "../envs/jupyter.yaml"
    shell:
        "jupyter nbconvert --to html {input}"


# rule convert_to_seurat:
#     input:
#         "{outdir}/preprocess/{stage}/{soloFeatures}.h5ad",
#     output:
#         "{outdir}/preprocess/{stage}/{soloFeatures}.rds",
#     conda:
#         "../envs/seurat.yaml"
