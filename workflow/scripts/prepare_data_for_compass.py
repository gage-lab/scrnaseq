#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = ["Michael Cuoco", "Joelle Faybishenko"]


import pegasus as pg
import sys
import pandas as pd
from pathlib import Path
from snakemake.shell import shell

sys.stderr = open(snakemake.log[0], "w")

print(f"opening file {snakemake.input[0]}")
data = pg.read_input(snakemake.input[0])

print("normalize data")
pg.normalize(data, norm_count=1e5)

outdir = Path(snakemake.output.norm).parent.parent
print(f"writing files to {outdir}")
pg.write_output(data, outdir, file_type="mtx")

for file in snakemake.output:
    print(f"decompressing {file}")
    shell(f"gunzip {file}.gz")

# fix barcodes and features files for proper reading by compass
df = pd.read_csv(snakemake.output.barcodes, sep="\t")
df["barcodekey"].to_csv(snakemake.output.barcodes, sep="\t", index=False, header=False)

df = pd.read_csv(snakemake.output.features, sep="\t")
df["featurekey"].to_csv(snakemake.output.features, sep="\t", index=False, header=False)


sys.stderr.close()
