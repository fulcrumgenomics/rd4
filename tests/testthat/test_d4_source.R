test_that("S3 methods don't break base/stats functions", {
  expect_equal(mean(c(1, 2, 3, 3)), 2.25)
  expect_equal(median(c(1, 2, 3, 3)), 2.5)
})

test_that("Methods raise error for nonexistent path", {
  # Note: helpful panic message from child thread is printed, but cannot be
  # captured here to test its contents.
  expect_error(seqinfo(D4Source("nonexistent_path")))
})

test_that("get_source() returns the original data path", {
  expect_equal(
    get_source(source_example1),
    testthat::test_path("testdata", "example1.d4")
  )
})

test_that("print() works for a D4Source object", {
  expect_equal(
    print(source_example1),
    paste(
      "<D4Source object>",
      paste("Source file:", testthat::test_path("testdata", "example1.d4")),
      "Number of chromosomes: 25",
      "Total chromosome length: 3088286401",
      "Number of tracks: 1",
      sep = "\n"
    )
  )
})

test_that(
  "seqinfo() returns the full sequence dictionary for files with one or
    multiple tracks",
  {
    exp_seq_dict <- Seqinfo(
      seqnames = c(
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM"
      ),
      seqlengths = c(
        248956422,
        242193529,
        198295559,
        190214555,
        181538259,
        170805979,
        159345973,
        145138636,
        138394717,
        133797422,
        135086622,
        133275309,
        114364328,
        107043718,
        101991189,
        90338345,
        83257441,
        80373285,
        58617616,
        64444167,
        46709983,
        50818468,
        156040895,
        57227415,
        16569
      )
    )
    expect_equal(seqinfo(source_example1), exp_seq_dict)
    expect_equal(seqinfo(source_multitrack), exp_seq_dict)
  }
)

test_that(
  "get_tracks() works for a D4 file with one track and no track label",
  {expect_equal(get_tracks(source_example1), c(""))}
)

test_that("get_tracks() works for a D4 file with multiple labeled tracks", {
  expect_equal(
    get_tracks(source_multitrack),
    c("track1_rep1", "track1_rep2", "track2")
  )
})

test_that("query() raises error for nonexistent track name", {
  # Note: helpful panic message from child thread is printed, but cannot be
  # captured here to test its contents.
  expect_error(query(source_example1, "chr7", 10, 20, "nonexistent_track"))
})

test_that("query() raises error for nonexistent contig name", {
  # Note: helpful panic message from child thread is printed, but cannot be
  # captured here to test its contents.
  expect_error(query(source_example1, "nonexistent_contig", 10, 20))
})

test_that("query() raises error for invalid coordinates", {
  # Note: helpful panic message from child thread is printed, but cannot be
  # captured here to test its contents.
  expect_error(query(source_example1, "chr7", 100, 98))
})

test_that("query() works for a region with no data", {
  expect_equal(query(source_example1, "chr7", 11, 20)$results, rep(0, 10))
})

test_that(
  "query() works for a region with some empty positions and some data",
  {
    query <- query(source_example1, "chr1", 43349437, 43349440)
    expect_equal(
      query$query,
      list(
        chr = "chr1",
        start = 43349437,
        end = 43349440,
        minus_strand = FALSE
      )
    )
    expect_equal(query$results, c(2, 1, 0, 0))
  }
)

test_that(
  "query() works for a minus strand region with some empty positions
    and some data",
  {
    query <-
      query(source_example1, "chr1", 43349437, 43349440,
            minus_strand = TRUE)
    expect_equal(query$query,
                 list(
                   chr = "chr1",
                   start = 43349437,
                   end = 43349440,
                   minus_strand = TRUE
                 ))
    expect_equal(query$results, c(0, 0, 1, 2))
  }
)

test_that("query() works for a region with multiple non-empty tracks", {
  expect_equal(
    query(source_multitrack, "chr1", 43349437, 43349440,
          track = "track1_rep2")$results,
    c(2, 1, 0, 0)
  )
})

test_that("to_granges() raises error for invalid chromosome", {
  expect_error(to_granges(source_multitrack, chr = "nonexistent_contig"))
})

test_that("to_granges() raises error for invalid coordinates", {
  expect_error(
    to_granges(source_multitrack, chr = "chrM", start = 100, end = 90))
})

test_that("to_granges() raises warning for coordinates with no chromosome", {
  expect_warning(to_granges(source_multitrack, start = 100, end = 200))
})

test_that("to_granges() raises warning for only one of start/end provided", {
  expect_warning(to_granges(source_multitrack, chr = "chrM", start = 100))
  expect_warning(to_granges(source_multitrack, chr = "chrM", end = 100))
})

test_that("to_granges() works for whole chromosome", {
  granges <- to_granges(source_multitrack, chr = "chr1")
  num_ranges <- length(granges)
  # Chromosome names
  expect_true(all(seqnames(granges) == rep("chr1", num_ranges)))
  # Interval coordinates
  expect_equal(ranges(granges)[1], IRanges(start = 1, end = 17027552))
  expect_equal(ranges(granges)[2], IRanges(start = 17027553, 17027562))
  expect_equal(
    ranges(granges)[num_ranges - 1], 
    IRanges(start = 193204137, end = 193204145))
  expect_equal(
    ranges(granges)[num_ranges], 
    IRanges(start = 193204146, end = 248956422))
  # Score values
  expect_equal(mcols(granges)$score[1], 0)
  expect_equal(mcols(granges)$score[2], 1)
  expect_equal(mcols(granges)$score[num_ranges - 1], 1)
  expect_equal(mcols(granges)$score[num_ranges], 0)
})

test_that("to_granges() works for specific region", {
  granges <- to_granges(
    source_multitrack, chr = "chr1", start = 17027561, end = 17027564)
  expect_equal(ranges(granges)[1], IRanges(start = 17027561, end = 17027562))
  expect_equal(ranges(granges)[2], IRanges(start = 17027563, end = 17027564))
  expect_equal(mcols(granges)$score[1], 1)
  expect_equal(mcols(granges)$score[2], 2)
})

test_that("to_granges() works for specific track", {
  granges <- to_granges(
    source_multitrack, 
    chr = "chr3", start = 37011633, 
    end = 37011640, 
    track = "track2"
  )
  expect_true(all(mcols(granges) == data.frame(score = 2)))
})

test_that("to_granges() works with custom metadata column name", {
  granges <- to_granges(source_multitrack, chr = "chrM", colname = "new_column")
  expect_true(all(mcols(granges) == data.frame(new_column = 0)))
})

test_that("import() works for whole genome", {
  granges <- import(source_multitrack)
  num_ranges <- length(granges)
  # Chromosome names
  expect_true(seqnames(granges)[1] == "chr1")
  expect_true(seqnames(granges)[num_ranges] == "chrM")
  # Interval coordinates
  expect_equal(ranges(granges)[1], IRanges(start = 1, end = 17027552))
  expect_equal(ranges(granges)[2], IRanges(start = 17027553, 17027562))
  expect_equal(ranges(granges)[num_ranges], IRanges(start = 1, end = 16569))
  # Score values
  expect_equal(mcols(granges)$score[1], 0)
  expect_equal(mcols(granges)$score[2], 1)
  expect_equal(mcols(granges)$score[num_ranges], 0)
})

test_that("import() works for one interval", {
  granges <- import(
    source_multitrack, 
    which = GRanges(
      seqnames = "chr1", 
      ranges = IRanges(start = 17027561, end = 17027564)))
  expect_equal(ranges(granges)[1], IRanges(start = 17027561, end = 17027562))
  expect_equal(ranges(granges)[2], IRanges(start = 17027563, end = 17027564))
  expect_equal(mcols(granges)$score[1], 1)
  expect_equal(mcols(granges)$score[2], 2)
})

test_that("import() simplifies overlapping query intervals", {
  granges <- import(
    source_multitrack, 
    which = GRanges(
      seqnames = "chr1", 
      ranges = IRanges(
        start = c(17027561, 17027562), 
        end = c(17027564, 17027563)))
    )
  expect_equal(ranges(granges)[1], IRanges(start = 17027561, end = 17027562))
  expect_equal(ranges(granges)[2], IRanges(start = 17027563, end = 17027564))
  expect_equal(mcols(granges)$score[1], 1)
  expect_equal(mcols(granges)$score[2], 2)
})

test_that("import() works for GRangesList with multiple chromosomes", {
  granges <- import(
    source_multitrack, 
    which = GRanges(
      seqnames = c("chr1", "chrM"), 
      ranges = IRanges(start = c(17027561, 1), end = c(17027564, 100))))
  expect_equal(as.vector(seqnames(granges)), c("chr1", "chr1", "chrM"))
  expect_equal(ranges(granges)[1], IRanges(start = 17027561, end = 17027562))
  expect_equal(ranges(granges)[2], IRanges(start = 17027563, end = 17027564))
  expect_equal(ranges(granges)[3], IRanges(start = 1, end = 100))
  expect_equal(mcols(granges)$score[1], 1)
  expect_equal(mcols(granges)$score[2], 2)
  expect_equal(mcols(granges)$score[3], 0)
})

test_that("import() works for specific track", {
  granges <- import(
    source_multitrack, 
    which = GRanges(
      seqnames = "chr3", 
      ranges = IRanges(start = 37011633, end = 37011640)),
    track = "track2"
  )
  expect_true(all(mcols(granges) == data.frame(score = 2)))
})

test_that("import() works with custom metadata column name", {
  granges <- import(
    source_multitrack, 
    which = GRanges(
      seqnames = "chrM", 
      ranges = IRanges(start = 1, end = 100)
    ), 
    colname = "new_column"
  )
  expect_true(all(mcols(granges) == data.frame(new_column = 0)))
})

test_that(
  "update_with_query_results() works for a GRanges with several test intervals",
  {
    # Test cases included:
    # - A region with no data
    # - A plus strand region with some empty positions and some data
    # - A region with some empty positions and some data and no strand
    # - A minus strand region with some empty positions and some data
    
    granges <-
      GenomicRanges::GRanges(
        seqnames = c("chr7", "chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(
          start = c(11, 43349437, 43349437, 43349437),
          end = c(20, 43349440, 43349440, 43349440)
        ),
        strand = c("+", "+", "*", "-")
      )
    granges_updated <-
      update_with_query_results(granges, source_example1)
    
    expect_true(identical(
      sapply(mcols(granges_updated)$query_results, function(x)
        x$results),
      list(rep(0, 10), c(2, 1, 0, 0), c(2, 1, 0, 0), c(0, 0, 1, 2))
    ))
    
    expect_true(identical(
      lapply(mcols(granges_updated)$query_results, function(x)
        x$query),
      list(
        list(
          chr = "chr7",
          start = 11,
          end = 20,
          minus_strand = FALSE
        ),
        list(
          chr = "chr1",
          start = 43349437,
          end = 43349440,
          minus_strand = FALSE
        ),
        list(
          chr = "chr1",
          start = 43349437,
          end = 43349440,
          minus_strand = FALSE
        ),
        list(
          chr = "chr1",
          start = 43349437,
          end = 43349440,
          minus_strand = TRUE
        )
      )
    ))
    
  })

test_that(
  "update_with_query_results() works for a region with multiple non-empty tracks",
  {
    granges <-
      GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
          start = 43349437, end = 43349440),
        strand = "+"
      )
    granges_updated <-
      update_with_query_results(granges, source_multitrack,
                                track = "track1_rep2")
    metadata <-
      mcols(granges_updated)$query_results_track1_rep2[[1]]
    expect_equal(metadata$results, c(2, 1, 0, 0))
    expect_equal(metadata$track, "track1_rep2")
  }
)

test_that(
  "mean() returns the same value as independently taking the mean of a
    query result",
  {
    expect_equal(
      mean(source_example2, chr = "chr3", start = 1, end = 5e+07),
      mean(query(source_example2, chr = "chr3", start = 1, end = 5e+07)$results)
    )
  }
)

test_that("mean() returns 0 for a region with no data", {
  expect_equal(
    mean(source_example2, chr = "chr3", start = 1e+06, end = 2e+06),
    0
  )
})

test_that("mean() works for a region with some empty positions and some data", {
  expect_equal(
    mean(source_example2, chr = "chr3", start = 37011631, end = 37011646), 
    2.125
  )
})

test_that("mean() works for a region with multiple non-empty tracks", {
  expect_equal(
    mean(
      source_multitrack, 
      chr = "chr3", 
      start = 37011716, 
      end = 37011720, 
      track = "track2"
    ),
    20.4
  )
})

test_that(
  "update_with_mean() works for a Granges with several test intervals", {
    granges <-
      GenomicRanges::GRanges(seqnames = "chr3", ranges = IRanges::IRanges(
        start = c(37011631, 1e+06),
        end = c(37011646, 2e+06)
      ))
    granges_updated <- update_with_mean(granges, source_example2)
    expect_true(identical(mcols(granges_updated)$mean, c(2.125, 0)))
  })

test_that(
  "update_with_mean() works with a user-specified metadata column name", {
    granges <-
      GenomicRanges::GRanges(seqnames = "chr3", ranges = IRanges::IRanges(
        start = c(37011631, 1e+06),
        end = c(37011646, 2e+06)
      ))
    granges_updated <- update_with_mean(
      granges, source_example2, colname = "newcol")
    expect_true(identical(mcols(granges_updated)$newcol, c(2.125, 0)))
  })

test_that(
  "update_with_mean() works for a region with multiple non-empty tracks",
  {
    granges <- GenomicRanges::GRanges(
      seqnames = "chr3",
      ranges = IRanges::IRanges(
        start = c(37011716), end = c(37011720))
    )
    granges_updated <-
      update_with_mean(granges, source_multitrack, track = "track2")
    expect_true(
      identical(mcols(granges_updated)$mean_track2, c(20.4))
    )
  }
)

test_that(
  "median() returns the same value as independently taking the median
    of a query result",
  {
    expect_equal(
      median(
        source_example2, chr = "chr3", start = 1, end = 5e+07),
      median(
        query(source_example2, chr = "chr3", start = 1,
              end = 5e+07)$results)
    )
  }
)

test_that("median() returns 0 for a region with no data", {
  expect_equal(
    median(source_multitrack, chr = "chr3", start = 1e+06, end = 2e+06,
           track = "track2"),
    0
  )
})

test_that(
  "median() works for a region with some empty positions and some data",
  {
    expect_equal(
      median(source_example2, chr = "chr3", start = 37011626,
             end = 37011640),
      2
    )
  }
)

test_that("median() works for a region with multiple non-empty tracks", {
  expect_equal(
    median(source_multitrack, chr = "chr3", start = 37011626,
           end = 37011640, track = "track2"),
    2
  )
})

test_that(
  "update_with_median() works for a Granges with several test intervals",
  {
    granges <-
      GenomicRanges::GRanges(seqnames = "chr3", ranges = IRanges::IRanges(
        start = c(1e+06, 37011626),
        end = c(2e+06, 37011640)
      ))
    granges_updated <-
      update_with_median(granges, source_example2)
    expect_true(identical(mcols(granges_updated)$median, c(0, 2)))
  }
)

test_that(
  "update_with_median() works with a user-specified column name",
  {
    granges <-
      GenomicRanges::GRanges(seqnames = "chr3", ranges = IRanges::IRanges(
        start = c(1e+06, 37011626),
        end = c(2e+06, 37011640)
      ))
    granges_updated <-
      update_with_median(granges, source_example2, colname = "newcol")
    expect_true(identical(mcols(granges_updated)$newcol, c(0, 2)))
  }
)

test_that(
  "update_with_median() works for a region with multiple non-empty tracks",
  {
    granges <-
      GenomicRanges::GRanges(
        seqnames = "chr3",
        ranges = IRanges::IRanges(
          start = c(37011626), end = c(37011640))
      )
    granges_updated <-
      update_with_median(granges, source_multitrack, track = "track2")
    expect_true(
      identical(mcols(granges_updated)$median_track2, c(2))
    )
  }
)

test_that("A histogram over a full range works", {
  hist <- source_example3$histogram("chr1", 0, 1000, NA, min = 0, max = NA)
  expect_equal(hist$mean(), 1.4)
  expect_equal(hist$median(), 0)
  expect_equal(hist$percentile(99.79), 100)
  expect_equal(hist$fraction_below(199), 0.998)
  expect_equal(hist$value_count(100), 2)
  expect_equal(hist$value_fraction(100), 0.002)
  expect_equal(hist$std(), 32.5275268042315)
})

test_that("A histogram over a partial range works", {
  hist <- source_example3$histogram("chr1", 0, 1000, NA, min = 99, max = 200)
  expect_equal(hist$mean(), 0.2)
  expect_equal(hist$median(), 99)
  expect_equal(hist$percentile(99.79), 100)
  expect_equal(hist$fraction_below(199), 0.998)
  expect_equal(hist$value_count(100), 2)
  expect_equal(hist$value_fraction(100), 0.002)
  expect_equal(hist$std(), 4.46766158073774)
})

test_that("percentile() works for a region with no data", {
  expect_equal(
    percentile(
      source_multitrack, "chr3", 1e+06, 2e+06, 100, track = "track2"),
    0
  )
})

test_that(
  "percentile() works for a region with some empty positions and some data",
  {
    expect_equal(
      percentile(source_multitrack, "chr3", 37011631, 37011646, 1,
                 track = "track2"),
      0
    )
    expect_equal(
      percentile(source_multitrack, "chr3", 37011631, 37011646, 50,
                 track = "track2"),
      2
    )
    expect_equal(
      percentile(source_multitrack, "chr3", 37011631, 37011646, 90,
                 track = "track2"),
      3
    )
    expect_equal(
      percentile(source_multitrack, "chr3", 37011631, 37011646, 100,
                 track = "track2"),
      3
    )
  }
)

test_that(
  "percentile() returns the same value as median() for percentile == 50",
  {
    expect_equal(
      percentile(source_multitrack, "chr3", 3.7e+07, 3.8e+07, 50,
                 track = "track2"),
      median(source_multitrack, chr = "chr3", start = 3.7e+07,
             end = 3.8e+07, track = "track2")
    )
  }
)

test_that("percentile() returns 0 for percentile <= 0", {
  expect_equal(
    percentile(source_multitrack, "chr3", 37011641, 37011646, 0,
               track = "track2"),
    0
  )
  expect_equal(
    percentile(source_multitrack, "chr3", 37011641, 37011646, -1,
               track = "track2"),
    0
  )
})

test_that("percentile() returns the max for percentile >= 100", {
  expect_equal(
    percentile(source_multitrack, "chr3", 37011631, 37011646, 100,
               track = "track2"),
    3
  )
  expect_equal(
    percentile(source_multitrack, "chr3", 37011631, 37011646, 101,
               track = "track2"),
    3
  )
})

test_that(
  "update_with_percentile() works for a Granges with several test intervals",
  {
    granges <-
      GenomicRanges::GRanges(
        seqnames = "chr3", ranges = IRanges::IRanges(
          start = c(1e+06, 37011626),
          end = c(2e+06, 37011640)
        ))
    granges_updated <-
      update_with_percentile(granges, source_example2, 50)
    expect_true(
      identical(mcols(granges_updated)$percentile_50, c(0, 2))
    )
  }
)

test_that(
  "update_with_percentile() works for a region with multiple non-empty tracks",
  {
    granges <- GenomicRanges::GRanges(
      seqnames = "chr3",
      ranges = IRanges::IRanges(
        start = c(37011626), end = c(37011640))
    )
    granges_updated <-
      update_with_percentile(
        granges,  
        source_multitrack,
        percentile = 50,
        track = "track2")
    expect_true(
      identical(
        mcols(granges_updated)$percentile_50_track2,
        c(2)
      )
    )
  }
)

test_that("Multiple calls to update Granges metadata behave as expected", {
  granges <-
    GenomicRanges::GRanges(seqnames = "chr3",
                           ranges = IRanges::IRanges(
                             start = c(37011626), end = c(37011640)))
  granges <-
    update_with_percentile(
      granges,  
      source_multitrack,
      percentile = 50,
      track = "track2")
  granges <-
    update_with_percentile(
      granges,  
      source_multitrack,
      percentile = 51,
      track = "track2")
  granges <-
    update_with_median(granges, source_multitrack, track = "track2")
  expect_true(identical(mcols(granges)$percentile_50_track2, c(2)))
  expect_true(identical(mcols(granges)$percentile_51_track2, c(2)))
  expect_true(identical(mcols(granges)$median_track2, c(2)))
})

test_that("resample() works for a region with no data", {
  expect_equal(
    resample(
      source_example1,
      "chr1",
      1e+06,
      2e+06,
      "median",
      bin_size = 1000,
      allow_bin_size_adjustment = FALSE
    )$results,
    rep(0, 1001)
  )
})

test_that(
  "resample() works for a region with some empty positions and some data and
    median method",
  {
    expect_equal(
      resample(
        source_example1,
        "chr1",
        17027541,
        17027570,
        "median",
        bin_size = 10,
        allow_bin_size_adjustment = FALSE
      )$results,
      c(0, 1, 2)
    )
  }
)

test_that("resample() works with default bin size (argument omitted)", {
  expect_equal(
    resample(source_example1, "chr1", 17027541, 17027570, "median")$results,
    c(1)
  )
})

test_that(
  "resample() works for a region with multiple non-empty tracks and mean
    method",
  {
    expect_equal(
      resample(
        source_multitrack,
        "chr1",
        17027541,
        17027570,
        "mean",
        bin_size = 10,
        allow_bin_size_adjustment = FALSE,
        track = "track1_rep1"
      )$results,
      c(0, 0.8, 2.2)
    )
  }
)
