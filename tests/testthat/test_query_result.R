result_example2 <- source_example2$query("chr3", 37011630, 37011640, NA)

test_that("results() returns the correct vector of values", {
  expect_equal(result_example2$results(), c(rep(0, 2), rep(2, 8)))
})

test_that("source() returns the correct data path", {
  expect_equal(result_example2$source(), testthat::test_path("testdata", "example2.d4"))
})

test_that("track() works when no track was provided", {
  expect_equal(result_example2$track(), "")
  expect_equal(source_multitrack$query("chr3", 37011630, 37011640, NA)$track(), "")
})

test_that("track() works when there are multiple tracks", {
  expect_equal(source_multitrack$query("chr3", 37011630, 37011640, "track2")$track(), "track2")
})

test_that("bin_size() returns NA if not resampled", {
  expect_true(is.na(result_example2$bin_size()))
})

test_that("bin_size() is set if resampled", {
  expect_equal(source_multitrack$resample("chr1", 17027540, 17027570, NA, "mean", 10, FALSE)$bin_size(), 10)
})
