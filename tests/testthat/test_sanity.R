test_that("1 + 1 is 2", {
    expect_equal(1 + 1, 2)
})

test_that("A D4 File can be opened and metadata extracted", {
    test_file <- testthat::test_path("testdata", "test.d4")
    file <- D4File$new(test_file)
    expect_equal(length(file$list_tracks()), 1)
    expect_equal(length(file$list_chroms()), 1)
    expect_equal(file$list_chroms(), c("chr1"))
})

test_that("A D4 File can be opened and queried", {
    test_file <- testthat::test_path("testdata", "test.d4")
    file <- D4File$new(test_file)
    result <- file$query("chr1", 12, 22)
    expect_equal(result$query()$chr(), "chr1")
    expect_equal(result$query()$left(), 12)
    expect_equal(result$query()$right(), 22)
    expect_equal(result$d4_file(), test_file)
    expect_equal(result$d4_track(), "")
    expect_equal(result$results(), c(0, 0, 0, 200, 0, 0, 0, 0, 100, 0))
})