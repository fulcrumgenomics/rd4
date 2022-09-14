# rd4 0.99.3

* Remove instructions to install from GitHub
* `print()` method for `D4Source` class
* Remove `data.table` dependency
* Replace `get_chroms()` with `seqinfo()` method for `D4Source`
* Replace 0-based half-open coordinates with 1-based closed coordinates; update tests
* Rename `update_*()` methods and move object being updated to first argument
* `update_with_mean()`, `update_with_median()`, and `update_with_percentile()` store the metadata in atomic columns instead of list columns
* `to_granges()` function to convert `D4Source` data to a `GRanges` object

# rd4 0.99.2

* Update INSTALL instructions for RPM-based Linux systems

# rd4 0.99.1

* Vignette conforms to Bioconductor guidelines
* Update system requirements in DESCRIPTION
* Add INSTALL
* Remove unused Python directory

# rd4 0.99.0

* Initial Bioconductor submission
