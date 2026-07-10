# mzmineOptimiser

`mzmineOptimiser` is a small R package that turns your MZmine parameter optimisation script into an installable workflow.

The package reads `config.toml`, runs the two-phase optimisation, writes the final optimised `.mzbatch` files, and saves CSV summaries and diagnostic plots.

## Install from GitHub

After pushing this folder to a GitHub repository, install it with:

```r
install.packages("remotes")
remotes::install_github("nava20ir/mzmineOptimiser")
```

## Create a config file

```r
library(mzmineOptimiser)
write_template_config("config.toml")
```

Edit `config.toml` for your project paths, QC files and parameter grid.

## Run

```r
library(mzmineOptimiser)
run_mzmine_optimisation("config.toml")
```

## TOML layout

The configuration file has four sections:

- `[paths]`: MZmine executable, base batch, output folders.
- `[qc]`: one QC file for Phase 1 and all QC files for Phase 2.
- `[optimization]`: number of random Phase 1 runs, number of top parameter sets for Phase 2.
- `[parameters]`: values used to build the full parameter grid.

## Notes

- This package assumes your base `.mzbatch` contains the MZmine module names used in the original script.
- It does not install MZmine. `mzmine_exe` in `config.toml` must point to your working MZmine executable.
- The package imports R dependencies through the `DESCRIPTION` file; it does not call `install.packages()` during runtime.
