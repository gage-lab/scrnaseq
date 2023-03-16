#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = ["Michael Cuoco", "Joelle Faybishenko"]

import numpy as np
import pandas as pd
import pegasus as pg
import tempfile, os
from pathlib import Path


# TODO: read data in with pegasus, normalize to counts per million, and save to a mtx file

print(f"opening file {snakemake.input[0]}")
data = pg.read_input(snakemake.input[0])
print("normalize data")
pg.normalize(data, norm_count=1e6)
print(f"writing file {snakemake.output[0]}")
pg.write_output(data, snakemake.output[0], file_type="mtx")
