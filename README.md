# scRNA-seq workflow for 10x data

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.16.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/gage-lab/scrnaseq/workflows/Tests/badge.svg?branch=main)](https://github.com/gage-lab/scrnaseq/actions?query=branch%3Amain+workflow%3ATests)

## Setup

### Testing

```bash
# generate 10x v3 scRNA-seq test data
cd .test/ngs-test-data
snakemake scrnaseq_10x_v3 --use-conda -c1 --conda-cleanup-pkgs cache

# test the workflow
cd ../..
snakemake all --use-conda -c1 --directory .test --show-failed-logs
```

### Configuration

1. Create two TSVs: (1) a "runsheet" that describes each 10x run and (2) a "patientsheet" that describes each patient

The runsheet should have the following columns

1. `run_id`: a unique identifier for each run
2. `patient_id`: a unique identifier for each patient
3. `r1`: the path to the R1 fastq file
4. `r2`: the path to the R2 fastq file

The patientsheet should have the following columns

1. `patient_id`: a unique identifier for each patient
2. `condition`: the condition of the patient

3. Enter the paths to the runsheet and patientsheet in the `config.yaml` file. Specify an output directory.

```yaml
outdir: <path to output directory>
runsheet: <path to runsheet>
patientsheet: <path to patientsheet>
```

### Running the workflow

```bash
# run the workflow all the way through
# -c16 = use 16 cores, this can be changed to whatever you want
snakemake all --use-conda -c16 --show-failed-logs

# run until a certain step
snakemake all --use-conda -c16 --show-failed-logs --until <rule name>

# IE: stop after the STARsolo
snakemake all --use-conda -c16 --show-failed-logs --until STARsolo

# IE: index the genome
snakemake STAR_index --use-conda -c16 --show-failed-logs
```
