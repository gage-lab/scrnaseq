#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

from pathlib import Path
from snakemake.shell import shell

# setup input fastqs and outdir
if type(snakemake.input.r1) is str:
    r1 = snakemake.input.r1
    r2 = snakemake.input.r2
else:
    assert len(snakemake.input.r1) == len(
        snakemake.input.r2
    ), "Unequal number of R1 and R2 files"
    r1 = ",".join(snakemake.input.r1)
    r2 = ",".join(snakemake.input.r2)

outdir = str(Path(snakemake.output.bam).parent) + "/"

# set readFilesIn and solo10xProtocol
print(f"10x chemistry: {snakemake.config['10x_chemistry']}")
if snakemake.config["10x_chemistry"] == "5prime":
    solo10xProtocol = "--soloBarcodeMate 1 --clip5pNbases 39 0 --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10"
    readFilesIn = f"--readFilesIn {r1} {r2}"
elif snakemake.config["10x_chemistry"] == "3prime_v2":
    solo10xProtocol = "--soloUMIlen 16 --clipAdapterType CellRanger4"
    readFilesIn = f"--readFilesIn {r2} {r1}"
elif snakemake.config["10x_chemistry"] == "3prime_v3":
    solo10xProtocol = "--soloUMIlen 12 --clipAdapterType CellRanger4"
    readFilesIn = f"--readFilesIn {r2} {r1}"
else:
    raise ValueError("Invalid 10x protocol")

# set read command
print(f"Read 1 file(s): {r1}")
print(f"Read 2 file(s): {r2}")
if r1.endswith(".gz"):
    readcmd = "gunzip -c"
elif r1.endswith(".bz2"):
    readcmd = "bunzip2 -c"
elif r1.endswith(".fastq") or r1.endswith(".fq"):
    readcmd = ""
else:
    raise ValueError("Invalid read file")

# set solobarcode read length
if snakemake.config["STARsolo"].get("soloBarcodeReadLength"):
    soloBarcodeReadLength = "--soloBarcodeReadLength " + str(
        snakemake.config["STARsolo"]["soloBarcodeReadLength"]
    )
else:
    soloBarcodeReadLength = ""

# set soloMultiMappers
if snakemake.config["STARsolo"].get("soloMultiMappers"):
    soloMultiMappers = "--soloMultiMappers " + str(
        snakemake.config["STARsolo"]["soloMultiMappers"]
    )
else:
    soloMultiMappers = ""

if snakemake.config["STARsolo"].get("outSAMmultNmax"):
    outSAMmultNmax = "--outSAMmultNmax " + str(
        snakemake.config["STARsolo"]["outSAMmultNmax"]
    )
else:
    outSAMmultNmax = ""

outFilterMultimapNmax = snakemake.config["STARsolo"]["outFilterMultimapNmax"]
winAnchorMultimapNmax = snakemake.config["STARsolo"]["winAnchorMultimapNmax"]

mem = str(50e9)  # use 50 GB of RAM for sorting

shell(
    "STAR "
    " --runThreadN {snakemake.threads}"
    " --genomeDir {snakemake.input.idx}"
    " --readFilesCommand {readcmd}"
    " {readFilesIn}"
    " --soloType CB_UMI_Simple"
    " {solo10xProtocol}"
    " {soloBarcodeReadLength}"
    " {soloMultiMappers}"
    " {outSAMmultNmax}"
    " --soloFeatures {snakemake.params.soloFeatures}"
    " --soloCellFilter EmptyDrops_CR"
    " --soloCBwhitelist {snakemake.input.whitelist}"
    " --soloOutFormatFeaturesGeneField3 -"
    " --soloOutFileNames outs genes.tsv barcodes.tsv matrix.mtx"
    " --outFilterMultimapNmax {outFilterMultimapNmax}"
    " --winAnchorMultimapNmax {winAnchorMultimapNmax}"
    " --limitBAMsortRAM {mem} "
    " --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM"
    " --outSAMtype BAM SortedByCoordinate"
    " --outFileNamePrefix {outdir}"
)
