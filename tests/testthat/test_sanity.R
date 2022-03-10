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
    result <- file$query("chr1", 12, 22, NA)
    expect_equal(result$query()$chr(), "chr1")
    expect_equal(result$query()$left(), 12)
    expect_equal(result$query()$right(), 22)
    expect_equal(result$source(), test_file)
    expect_equal(result$track(), "")
    expect_equal(result$results(), c(0, 0, 0, 200, 0, 0, 0, 0, 100, 0))
})