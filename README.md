# rd4

## Overview

The Dense Depth Data Dump (D4) format provides fast analysis and compact storage of quantitative genomics datasets. `rd4` provides R bindings for reading and querying D4 files. For full details on the format, see Hou et al. (https://doi.org/10.1038/s43588-021-00085-0).

## Installation

TBA

## Usage

```{r}
library(rd4)
browseVignettes("rd4")
```

## For developers

### R dependencies

Consider using [renv](https://rstudio.github.io/renv/articles/renv.html) to isolate an R environment for `rd4` development.

- [remotes](https://cran.r-project.org/web/packages/remotes/index.html)
- [rextendr](https://github.com/extendr/rextendr)
- [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html)
- [BiocCheck](https://bioconductor.org/packages/release/bioc/html/BiocCheck.html)

### Developing an R package using Rust code

Note: these steps were already done for `rd4` and don't need to be repeated by new developers.

- [rextendr](https://extendr.github.io/rextendr/index.html) package
- [Instructions](https://extendr.github.io/rextendr/articles/package.html) on setting up and developing a package


### Build and test the package locally

Build Rust code

```bash
cd src/rust/
cargo build
```

Build and check R package within an R session

```R
# Set working directory to package root
setwd(".")

# Compile Rust code into R functions and auto-generate R documentation (yes, rextendr::document() does both)
rextendr::document()

# Load the package
devtools::load_all()

# Run tests
devtools::test()

# Run CRAN and/or Bioconductor checks
devtools::check()
BiocCheck::BiocCheck()
```


### Add new Rust code

[rextendr vignette](https://extendr.github.io/rextendr/articles/package.html)


### Test Rust code

```bash
cd src/rust/
cargo test
```

### Formatting of auto-generated R code

Auto-generated R code produced by `rextendr` may not satisfy formatting requirements for Bioconductor or CRAN. Consider using an automatic formatter like [formatr](https://yihui.org/formatr/) to auto-format before submission, or the `Code -> Reformat Code` workflow in Rstudio.
