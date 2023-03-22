#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = ["Michael Cuoco", "Joelle Faybishenko"]


import pegasus as pg
import sys
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

sys.stderr.close()
