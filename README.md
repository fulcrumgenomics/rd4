# rd4

## Overview

TBA

## Installation

TBA

## Usage

```R
# Open a D4 file
d4_file <- "file.d4"
d4_source <- D4Source$new(d4_file)

# Identify chromosomes and tracks
chroms <- d4_source$get_chroms()
tracks <- d4_source$get_tracks()

# Query a region
query <- d4_source$query("chr1", 100, 100000, NA)
query_results <- query$results()
query_chr <- query$query()$chr()

# Summary statistics for a region
mean <- d4_source$mean("chr1", 100, 100000, NA)
median <- d4_source$median("chr1", 100, 100000, NA)
histogram <- d4_source$histogram("chr1", 100, 100000, "", min = 0, max = 10)
histogram_total <- histogram$total_count()
percentile <- d4_source$percentile("chr1", 100, 100000, "", percentile = 0.5)

# Resample a region
resample <- d4_source$resample("chr1", 50000, 100000, NA, method = "mean", bin_size = 10, allow_bin_size_adjustment = FALSE)
resample_results <- resample$results()
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

```R
# Set working directory to package root
setwd(".")

# Compile Rust code into R functions and auto-generate R documentation (yes, rextendr::document() does both)
rextendr::document()

# Load the package
devtools::load_all(".")

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
