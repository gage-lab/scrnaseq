# CONTRIBUTE

## Setup your environment

```bash
# clone the repo
git clone --recurse-submodules https://github.com/gage-lab/scrnaseq.git

# create a development environment
mamba env create -f environment.yaml

# install the pre-commit hooks
pre-commit install
```

## Develop a new feature

1. See issues for a list of things to do.
2. Choose an issue to work on
3. Create a branch for the issue
4. Make changes and commit them
5. Create a pull request to activate the testing workflow in Github Actions
6. Make additional changes if necessary, until tests pass
7. Once tests pass, request code review from @mikecuoco
8. Once reviews are all addressed, @mikecuoco will merge the pull request to the main branch!

## Notebook development tips

1. When creating a rule with a notebook, edit the json of the notebook so that `"nbformat_minor": 5`
2. You can edit the run notebook in it's output location. **Be sure to select the conda environment that was created by snakemake as the jupyter kernel.**
3. Now you can develop in that notebook and copy-paste your code to the original notebook in the `workflow/notebooks` directory.
