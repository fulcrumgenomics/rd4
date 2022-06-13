Sys.setenv(TESTTHAT_MAX_FAILS = Inf)
source_example1 <- D4Source(testthat::test_path("testdata", "example1.d4"))
source_example2 <- D4Source(testthat::test_path("testdata", "example2.d4"))
source_example3 <- D4Source(testthat::test_path("testdata", "example3.d4"))
source_multitrack <- D4Source(
    testthat::test_path("testdata", "example_multitrack.d4"))
