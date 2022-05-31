test_that("get_source() returns the original data path", {
  expect_equal(source_example1$get_source(), testthat::test_path("testdata", "example1.d4"))
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
  expect_equal(source_example1$get_chroms(), exp_seq_dict)
  expect_equal(source_multitrack$get_chroms(), exp_seq_dict)
})

test_that("get_tracks() works for a D4 file with one track and no track label", {
  expect_equal(source_example1$get_tracks(), c(""))
})

test_that("get_tracks() works for a D4 file with multiple labeled tracks", {
  expect_equal(source_multitrack$get_tracks(), c("track1_rep1", "track1_rep2", "track2"))
})

test_that("query() raises error for nonexistent track name", {
  expect_error(source_example1$query("chr7", 10, 20, "nonexistent_track"), "Track does not exist: nonexistent_track")
})

test_that("query() raises error for nonexistent contig name", {
  expect_error(source_example1$query("nonexistent_contig", 10, 20, NA), "Contig does not exist: nonexistent_contig")
})

test_that("query() raises error for invalid coordinates", {
  expect_error(source_example1$query("chr1", -1, 10, NA), "Invalid coordinates: -1-10")
  expect_error(source_example1$query("chr7", 100, 999999999, NA), "Invalid coordinates: 100-999999999")
  expect_error(source_example1$query("chr7", 100, 99, NA), "Invalid coordinates: 100-99")
})

test_that("query() raises error for missing track name when there are multiple tracks", {
  expect_error(
    source_multitrack$query("chr1", 43349436, 43349440, NA), 
    "Track name required for D4 file with multiple tracks"
  )
})

test_that("query() works for a region with no data", {
  expect_equal(source_example1$query("chr7", 10, 20, NA)$results(), rep(0, 10))
})

test_that("query() returns an empty vector for a region of length 0", {
  expect_equal(source_example1$query("chr7", 10, 10, NA)$results(), c())
})

test_that("query() works for a region with some empty positions and some data", {
  expect_equal(source_example1$query("chr1", 43349436, 43349440, NA)$results(), c(2, 1, 0, 0))
})

test_that("query() works for a region with multiple non-empty tracks", {
  expect_equal(source_multitrack$query("chr1", 43349436, 43349440, "track1_rep2")$results(), c(2, 1, 0, 0))
})

test_that("mean() returns the same value as independently taking the mean of a query result", {
  expect_equal(
    source_example2$mean("chr3", 37000000, 38000000, NA),
    mean(source_example2$query("chr3", 37000000, 38000000, NA)$results()),
  )
})

test_that("mean() returns 0 for a region with no data", {
  expect_equal(source_example2$mean("chr3", 1000000, 2000000, NA), 0)
})

test_that("mean() works for a region with some empty positions and some data", {
  expect_equal(source_example2$mean("chr3", 37011630, 37011646, NA), 2.125)
})

test_that("mean() works for a region with multiple non-empty tracks", {
  expect_equal(source_multitrack$mean("chr3", 37011715, 37011720, "track2"), 20.4)
})

test_that("median() returns the same value as independently taking the median of a query result", {
  expect_equal(
    source_example2$median("chr3", 37000000, 38000000, NA),
    median(source_example2$query("chr3", 37000000, 38000000, NA)$results()),
  )
})

test_that("median() returns 0 for a region with no data", {
  expect_equal(source_multitrack$mean("chr3", 1000000, 2000000, "track2"), 0)
})

test_that("median() works for a region with some empty positions and some data", {
  expect_equal(source_example2$median("chr3", 37011624, 37011640, NA), 1)
})

test_that("median() works for a region with multiple non-empty tracks", {
  expect_equal(source_multitrack$median("chr3", 37011624, 37011640, "track2"), 1)
})

# test_that("histogram() raises error for invalid bucket limits", {expect_true(FALSE)})
# 
# test_that("histogram() works with default min and max bucket", {expect_true(FALSE)})
# 
# test_that("histogram() returns entire chromosome when coordinates are omitted: chromosome with data", {expect_true(FALSE)})
# 
# test_that("histogram() returns entire chromosome when coordinates are omitted: chromosome with no data", {expect_true(FALSE)})
# 
# test_that("histogram() works for a region with no data", {expect_true(FALSE)})
# 
# test_that("histogram() works for a region with some empty positions and some data", {expect_true(FALSE)})
# 
# test_that("histogram() works for a region with multiple non-empty tracks", {expect_true(FALSE)})

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
  hist <- source_example3$histogram("chr1", 0, 1000, NA, min=99, max=200)
  expect_equal(hist$mean(), 0.2)
  expect_equal(hist$median(), 99)
  expect_equal(hist$percentile(99.79), 100)
  expect_equal(hist$fraction_below(199), 0.998)
  expect_equal(hist$value_count(100), 2)
  expect_equal(hist$value_fraction(100), 0.002)
  expect_equal(hist$std(), 4.467661580737736)
})

test_that("percentile() raises error for invalid percentile value", {expect_true(FALSE)})

test_that("percentile() works for a region with no data", {
  expect_equal(source_multitrack$percentile("chr3", 1000000, 2000000, "track2", 1), 0)
})

test_that("percentile() works for a region with some empty positions and some data", {
  expect_equal(source_multitrack$percentile("chr3", 37011630, 37011646, "track2", 0.1), 0)
  expect_equal(source_multitrack$percentile("chr3", 37011630, 37011646, "track2", 0.5), 2)
  expect_equal(source_multitrack$percentile("chr3", 37011630, 37011646, "track2", 0.9), 3)
  expect_equal(source_multitrack$percentile("chr3", 37011630, 37011646, "track2", 1), 3)
})

test_that("percentile() returns the same value as median() for percentile=0.5", {
  expect_equal(
    source_multitrack$percentile("chr3", 37000000, 38000000, "track2", 0.5), 
    source_multitrack$median("chr3", 37000000, 38000000, "track2")
  )
})

test_that("A D4 File can be opened and a region can be resampled", {
  result <- source_example3$resample(
    "chr1", 12, 22, track=NA, method="mean", bin_size=NA, 
    allow_bin_size_adjustment=NA)
  expect_equal(length(result$results()), 1)
  
  result <- source_example3$resample(
    "chr1", 12, 22, track=NA, method="median", bin_size=NA, 
    allow_bin_size_adjustment=NA)
  expect_equal(length(result$results()), 1)
  
  result <- source_example3$resample(
    "chr1", 0, 1000, track=NA, method="mean", bin_size=10, 
    allow_bin_size_adjustment=NA)
  expect_equal(length(result$results()), 100)
  
  result <- source_example3$resample(
    "chr1", 0, 1000, track=NA, method="median", bin_size=10, 
    allow_bin_size_adjustment=NA)
  expect_equal(length(result$results()), 100)
})

test_that("Methods with a track parameter work with default missing value", {
  expect_equal(source_example1$query("chr1", 43349436, 43349440)$results(), c(2, 1, 0, 0))
  expect_equal(source_multitrack$query("chr1", 43349436, 43349440)$results(), c(2, 1, 0, 0))
  expect_equal(source_example2$mean("chr3", 37011630, 37011646), 2.125)
  expect_true(FALSE) # mean() multitrack
  expect_true(FALSE) # median() single track
  expect_true(FALSE) # median() multitrack
  # expect_true(FALSE) # histogram() single track
  # expect_true(FALSE) # histogram() multitrack
  expect_true(FALSE) # percentile() single track
  expect_true(FALSE) # percentile() multitrack
  expect_true(FALSE) # resample() single track
  expect_true(FALSE) # resample() multitrack
})
