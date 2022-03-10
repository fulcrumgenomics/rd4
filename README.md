# rd4

## Install R Deps

- R itself
    - https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-20-04-quickstart
- [`usethis`](https://usethis.r-lib.org/)
    - `install.packages("usethis")`
- [`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
    - `install.packages("devtools")`
- [`remotes`]()
    - `install.packages("remotes")`
- [`rextendr`]()
    - `remotes::install_github("extendr/rextendr")`

## Create the template repo

```R
usethis::create_package("~/dev/rd4")
usethis::use_testthat()
setwd("~/dev/rd4")
rextendr::use_extendr()
rextendr::document()
```

## Build and test things manually

```R
setwd("~/dev/rd4")
rextendr::document()
devtools::load_all(".")
```

### Test Rust code

```bash
cd src/rust/
cargo test
```

### Test R Code

```R
# In development dir
devtools::check() # Will check that the package as a whole is well formed
devtools::test() # will run tets in `tests/testthat`
```

## Example

```R
devtools::load_all(".")
file <- D4File::new("path_to_file.d4")
result <- file$query("chr1", 100, 1000)
result <- result$result()
```

