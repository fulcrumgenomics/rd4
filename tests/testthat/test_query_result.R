result_example2 <- query(source_example2, "chr3", 37011630, 37011640)

test_that("results() returns the correct vector of values", {
  expect_equal(results(result_example2), c(rep(0, 2), rep(2, 8)))
})

test_that("source() returns the correct data path", {
  expect_equal(source(result_example2), testthat::test_path("testdata", "example2.d4"))
})

test_that("track() works when no track was provided", {
  expect_equal(track(result_example2), "")
  expect_equal(track(query(source_multitrack, "chr3", 37011630, 37011640, NA)), "")
})

test_that("track() works when there are multiple tracks", {
  expect_equal(track(query(source_multitrack, "chr3", 37011630, 37011640, "track2")), "track2")
})

test_that("bin_size() returns NA if not resampled", {
  expect_true(is.na(bin_size(result_example2)))
})

test_that("bin_size() is set if resampled", {
  expect_equal(bin_size(resample(source_multitrack, "chr1", 17027540, 17027570, NA, "mean", 10, FALSE)), 10)
})
