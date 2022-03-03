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
setwd("~/dev/rd4")
rextendr::use_extendr()
rextendr::document()
```

## Build and test things

```R
setwd("~/dev/rd4")
rextendr::document()
devtools::load_all(".")
```
