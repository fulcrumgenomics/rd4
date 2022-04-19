test_that("A D4 File can be opened and metadata extracted", {
    test_file <- testthat::test_path("testdata", "test.d4")
    file <- D4Source$new(test_file)
    expect_equal(length(file$get_tracks()), 1)
    expect_equal(length(file$get_chroms()), 1)
    expect_equal(file$get_chroms(), list(list("name"="chr1", "size"=1000)))
})

test_that("A D4 File can be opened and queried", {
    test_file <- testthat::test_path("testdata", "test.d4")
    file <- D4Source$new(test_file)

    expect_equal(file$mean("chr1", 12, 22, track=NA), 30.0)
    expect_equal(file$median("chr1", 12, 22, track=NA), 0.0)
    expect_equal(file$percentile("chr1", 12, 22, track=NA, 50.0), 0.0)

    result <- file$query("chr1", 12, 22, NA)
    expect_equal(result$query()$chr(), "chr1")
    expect_equal(result$query()$left(), 12)
    expect_equal(result$query()$right(), 22)
    expect_equal(result$source(), test_file)
    expect_equal(result$track(), "")
    expect_equal(result$results(), c(0, 0, 0, 200, 0, 0, 0, 0, 100, 0))
})

test_that("A D4 File can be opened and a region can be resampled", {
    test_file <- testthat::test_path("testdata", "test.d4")
    file <- D4Source$new(test_file)

    result <- file$resample("chr1", 12, 22, track=NA, method="mean", bin_size=NA, allow_bin_size_adjustment=NA)
    expect_equal(length(result$results()), 1)

    result <- file$resample("chr1", 12, 22, track=NA, method="median", bin_size=NA, allow_bin_size_adjustment=NA)
    expect_equal(length(result$results()), 1)

    result <- file$resample("chr1", 0, 1000, track=NA, method="mean", bin_size=10, allow_bin_size_adjustment=NA)
    expect_equal(length(result$results()), 100)

    result <- file$resample("chr1", 0, 1000, track=NA, method="median", bin_size=10, allow_bin_size_adjustment=NA)
    expect_equal(length(result$results()), 100)
})

test_that("A Histogram over a full range works", {
    test_file <- testthat::test_path("testdata", "test.d4")
    file <- D4Source$new(test_file)
    hist <- file$histogram("chr1", 0, 1000, NA, min=0, max=NA)
    expect_equal(hist$mean(), 1.4)
    expect_equal(hist$median(), 0.0)
    expect_equal(hist$percentile(99.79), 100)
    expect_equal(hist$fraction_below(199), 0.998)
    expect_equal(hist$value_count(100), 2)
    expect_equal(hist$value_fraction(100), 0.002)
    expect_equal(hist$std(), 32.52752680423152)
})

test_that("A Histogram over a partial range works", {
    test_file <- testthat::test_path("testdata", "test.d4")
    file <- D4Source$new(test_file)

    hist <- file$histogram("chr1", 0, 1000, NA, min=99, max=200)
    expect_equal(hist$mean(), 0.2)
    expect_equal(hist$median(), 99)
    expect_equal(hist$percentile(99.79), 100)
    expect_equal(hist$fraction_below(199), 0.998)
    expect_equal(hist$value_count(100), 2)
    expect_equal(hist$value_fraction(100), 0.002)
    expect_equal(hist$std(), 4.467661580737736)
})

