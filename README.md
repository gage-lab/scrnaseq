# scRNA-seq workflow for 10x data

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.16.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/gage-lab/scrnaseq/workflows/Tests/badge.svg?branch=main)](https://github.com/gage-lab/scrnaseq/actions?query=branch%3Amain+workflow%3ATests)

## Usage

```bash
# generate 10x v3 scRNA-seq test data
snakemake scrnaseq_10x_v3 --use-conda -c1 --conda-cleanup-pkgs cache

# test the workflow
snakemake all --use-conda -c1 --directory .test
```

## Development
