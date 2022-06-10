result_example2 <- query(source_example2, "chr3", 37011630, 37011640)

test_that("results is the correct vector of values", {
  expect_equal(result_example2$results, c(rep(0, 2), rep(2, 8)))
})

test_that("source is the correct data path", {
  expect_equal(result_example2$source, testthat::test_path("testdata", "example2.d4"))
})

test_that("track works when no track was provided", {
  expect_equal(result_example2$track, "")
  expect_equal(query(source_multitrack, "chr3", 37011630, 37011640, track = NA)$track, "")
})

test_that("track works when there are multiple tracks", {
  expect_equal(query(source_multitrack, "chr3", 37011630, 37011640, track = "track2")$track, "track2")
})

test_that("bin_size is NA if not resampled", {
  expect_true(is.na(result_example2$bin_size))
})

test_that("bin_size is set if resampled", {
  expect_equal(resample(
    source_multitrack, "chr1", 17027540, 17027570, "mean", bin_size = 10, allow_bin_size_adjustment = FALSE)$bin_size, 
    10)
})
