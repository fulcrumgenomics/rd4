test_that("Methods raise error for nonexistent path", {
  # Note: helpful panic message from child thread is printed, but cannot be captured here to test its contents.
  expect_error(get_chroms(D4Source("nonexistent_path")))
})

test_that("get_source() returns the original data path", {
  expect_equal(get_source(source_example1), testthat::test_path("testdata", "example1.d4"))
})

test_that("get_chroms() returns the full sequence dictionary for files with one or multiple tracks", {
  exp_seq_dict <- list(
    list(name="chr1", size=248956422),
    list(name="chr2", size=242193529),
    list(name="chr3", size=198295559),
    list(name="chr4", size=190214555),
    list(name="chr5", size=181538259),
    list(name="chr6", size=170805979),
    list(name="chr7", size=159345973),
    list(name="chr8", size=145138636),
    list(name="chr9", size=138394717),
    list(name="chr10", size=133797422),
    list(name="chr11", size=135086622),
    list(name="chr12", size=133275309),
    list(name="chr13", size=114364328),
    list(name="chr14", size=107043718),
    list(name="chr15", size=101991189),
    list(name="chr16", size=90338345),
    list(name="chr17", size=83257441),
    list(name="chr18", size=80373285),
    list(name="chr19", size=58617616),
    list(name="chr20", size=64444167),
    list(name="chr21", size=46709983),
    list(name="chr22", size=50818468),
    list(name="chrX", size=156040895),
    list(name="chrY", size=57227415),
    list(name="chrM", size=16569)
  )
  expect_equal(get_chroms(source_example1), exp_seq_dict)
  expect_equal(get_chroms(source_multitrack), exp_seq_dict)
})

test_that("get_tracks() works for a D4 file with one track and no track label", {
  expect_equal(get_tracks(source_example1), c(""))
})

test_that("get_tracks() works for a D4 file with multiple labeled tracks", {
  expect_equal(get_tracks(source_multitrack), c("track1_rep1", "track1_rep2", "track2"))
})

test_that("query() raises error for nonexistent track name", {
  # Note: helpful panic message from child thread is printed, but cannot be captured here to test its contents.
  expect_error(query(source_example1, "chr7", 10, 20, "nonexistent_track"))
})

test_that("query() raises error for nonexistent contig name", {
  # Note: helpful panic message from child thread is printed, but cannot be captured here to test its contents.
  expect_error(query(source_example1, "nonexistent_contig", 10, 20))
})

test_that("query() raises error for invalid coordinates", {
  # Note: helpful panic message from child thread is printed, but cannot be captured here to test its contents.
  expect_error(query(source_example1, "chr7", 100, 99))
})

test_that("query() returns an empty result for a 0-length region", {
  expect_equal(query(source_example1, "chr7", 10, 10)$results, numeric())
})

test_that("query() works for a region with no data", {
  expect_equal(query(source_example1, "chr7", 10, 20)$results, rep(0, 10))
})

test_that("query() works for a region with some empty positions and some data", {
  expect_equal(query(source_example1, "chr1", 43349436, 43349440)$results, c(2, 1, 0, 0))
})

test_that("query() works for a region with multiple non-empty tracks", {
  expect_equal(query(source_multitrack, "chr1", 43349436, 43349440, "track1_rep2")$results, c(2, 1, 0, 0))
})

test_that("mean() returns the same value as independently taking the mean of a query result", {
  expect_equal(
    mean(source_example2, "chr3", 0, 50000000),
    mean(query(source_example2, "chr3", 0, 50000000)$results),
  )
})

test_that("mean() returns 0 for a region with no data", {
  expect_equal(mean(source_example2, "chr3", 1000000, 2000000), 0)
})

test_that("mean() works for a region with some empty positions and some data", {
  expect_equal(mean(source_example2, "chr3", 37011630, 37011646), 2.125)
})

test_that("mean() works for a region with multiple non-empty tracks", {
  expect_equal(mean(source_multitrack, "chr3", 37011715, 37011720, "track2"), 20.4)
})

test_that("median() returns the same value as independently taking the median of a query result", {
  expect_equal(
    median(source_example2, "chr3", 0, 50000000),
    median(query(source_example2, "chr3", 0, 50000000)$results)
  )
})

test_that("median() returns 0 for a region with no data", {
  expect_equal(median(source_multitrack, "chr3", 1000000, 2000000, "track2"), 0)
})

test_that("median() works for a region with some empty positions and some data", {
  expect_equal(median(source_example2, "chr3", 37011625, 37011640), 2)
})

test_that("median() works for a region with multiple non-empty tracks", {
  expect_equal(median(source_multitrack, "chr3", 37011625, 37011640, "track2"), 2)
})

test_that("A histogram over a full range works", {
  hist <- source_example3$histogram("chr1", 0, 1000, NA, min=0, max=NA)
  expect_equal(hist$mean(), 1.4)
  expect_equal(hist$median(), 0.0)
  expect_equal(hist$percentile(99.79), 100)
  expect_equal(hist$fraction_below(199), 0.998)
  expect_equal(hist$value_count(100), 2)
  expect_equal(hist$value_fraction(100), 0.002)
  expect_equal(hist$std(), 32.52752680423152)
})

test_that("A histogram over a partial range works", {
  hist <- source_example3$histogram("chr1", 0, 1000, NA, min=99, max=200)
  expect_equal(hist$mean(), 0.2)
  expect_equal(hist$median(), 99)
  expect_equal(hist$percentile(99.79), 100)
  expect_equal(hist$fraction_below(199), 0.998)
  expect_equal(hist$value_count(100), 2)
  expect_equal(hist$value_fraction(100), 0.002)
  expect_equal(hist$std(), 4.467661580737736)
})

test_that("percentile() works for a region with no data", {
  expect_equal(percentile(source_multitrack, "chr3", 1000000, 2000000, 100, "track2"), 0)
})

test_that("percentile() works for a region with some empty positions and some data", {
  expect_equal(percentile(source_multitrack, "chr3", 37011630, 37011646, 1, "track2"), 0)
  expect_equal(percentile(source_multitrack, "chr3", 37011630, 37011646, 50, "track2"), 2)
  expect_equal(percentile(source_multitrack, "chr3", 37011630, 37011646, 90, "track2"), 3)
  expect_equal(percentile(source_multitrack, "chr3", 37011630, 37011646, 100, "track2"), 3)
})

test_that("percentile() returns the same value as median() for percentile == 50", {
  expect_equal(
    percentile(source_multitrack, "chr3", 37000000, 38000000, 50, "track2"), 
    median(source_multitrack, "chr3", 37000000, 38000000, "track2")
  )
})

test_that("percentile() returns 0 for percentile <= 0", {
  expect_equal(percentile(source_multitrack, "chr3", 37011640, 37011646, 0, "track2"), 0)
  expect_equal(percentile(source_multitrack, "chr3", 37011640, 37011646, -1, "track2"), 0)
})

test_that("percentile() returns the max for percentile >= 100", {
  expect_equal(percentile(source_multitrack, "chr3", 37011630, 37011646, 100, "track2"), 3)
  expect_equal(percentile(source_multitrack, "chr3", 37011630, 37011646, 101, "track2"), 3)
})

test_that("resample() works for a region with no data", {
  expect_equal(resample(source_example1, "chr1", 1000000, 2000000, "median", 1000, FALSE)$results, rep(0, 1000))
})

test_that("resample() works for a region with some empty positions and some data and median method", {
  expect_equal(resample(source_example1, "chr1", 17027540, 17027570, "median", 10, FALSE)$results, c(0, 1, 2))
})

test_that("resample() works with default bin size (argument omitted)", {
  expect_equal(resample(source_example1, "chr1", 17027540, 17027570, "median")$results, c(1))
})

test_that("resample() works for a region with multiple non-empty tracks and mean method", {
  expect_equal(
    resample(source_multitrack, "chr1", 17027540, 17027570, "mean", 10, FALSE, "track1_rep1")$results, 
    c(0, 0.8, 2.2)
  )
})
