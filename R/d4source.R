# Convert a Rust QueryResultEnv to an R list
query_result <- function(query_result_env, minus_strand = FALSE) {
    results <- query_result_env$results()
    list(
        results = if (minus_strand)
            rev(results)
        else
            results,
        query = list(
            chr = query_result_env$query()$chr(),
            start = query_result_env$query()$left(),
            end = query_result_env$query()$right(),
            minus_strand = minus_strand
        ),
        source = query_result_env$source(),
        track = query_result_env$track(),
        bin_size = query_result_env$bin_size()
    )
}

#' Establish a connection to a D4 file
#'
#' @param path File path
#' @return Object of class `D4Source`
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' methods(class = 'D4Source')
#'
#' @export
D4Source <- function(path) {
    d4 <- D4SourceEnv$new(path)
    class(d4) <- c("D4Source", class(d4))
    d4
}

#' Returns the path to the source D4 file
#'
#' @param d4source `D4Source` object
#' @return The path to the source D4 file
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' get_source(d4source)
#'
#' @export
get_source <- function(d4source) {
    d4source$get_source()
}

#' Returns the list of chromosomes in the source D4 file
#'
#' @param d4source `D4Source` object
#' @return Nested list; each element represents a chromosome and stores the
#' chromosome name and length
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' get_chroms(d4source)
#'
#' @export
get_chroms <- function(d4source) {
    d4source$get_chroms()
}

#' Returns the names of tracks in the source D4 file
#'
#' @param d4source `D4Source` object
#' @return Vector of track names
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' get_tracks(d4source)
#'
#' @export
get_tracks <- function(d4source) {
    d4source$get_tracks()
}

#' Extracts the given region from the D4 file
#'
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param ... Not used for values, forces later arguments to bind by name
#' @param minus_strand If TRUE, vector of position scores will read in the
#' minus strand direction from
#' `end - 1` to `start`.
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Query result as a list
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' query_result <- query(d4source, 'chr1', 1000, 2000)
#'
#' @importFrom wrapr stop_if_dot_args
#' @export
query <-
    function(d4source, chr, start, end, ..., minus_strand = FALSE, track = NA) {
        stop_if_dot_args(substitute(list(...)), "query()")
        query_result(
            d4source$query(chr, start, end, track), minus_strand = minus_strand)
    }

#' For a set of genomic ranges in a `GRanges` object, extracts each region from
#' the D4 file and adds the query result as a metadata attribute for each
#' genomic range, returning a new `GRanges` object with updated metadata.
#'
#' @param d4source `D4Source` object
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @return New version of `granges` with query results added as the
#' `'query_results'` metadata attribute for each genomic range. For ranges on
#' the minus strand, the vector of position scores reads in the minus strand
#' direction. If a track is specified, the track name will be included in the
#' metadata column name.
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' granges <- GRanges(seqnames = c('chr1'), ranges = IRanges(1000:2000), 
#'                    strand = c('+'))
#' granges_with_query <- update_query_results(d4source, granges)
#'
#' @importFrom GenomicRanges mcols seqnames start end ranges strand
#' @importFrom wrapr stop_if_dot_args
#' @export
update_query_results <-
    function(d4source, granges, ..., track = NA) {
        stop_if_dot_args(substitute(list(...)), "update_query_results()")
        fn <- function(chr, start, end, minus_strand, track) {
            query(
                d4source,
                chr = chr,
                start = start,
                end = end,
                minus_strand = minus_strand,
                track = track
            )
        }
        mcol_name <- if (is.na(track))
            "query_results"
        else
            paste("query_results", track, sep = "_")
        GenomicRanges::mcols(granges)[[mcol_name]] <-
            mapply(
                fn,
                chr = as.vector(seqnames(granges)),
                start = as.vector(start(ranges(granges))) - 1,
                end = as.vector(end(ranges(granges))),
                minus_strand = as.vector(strand(granges) == "-"),
                track = track,
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )

        granges

    }

#' Returns the mean data value over the given region
#'
#' @param x `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Mean data value
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' mean(d4source, 'chr1', 1000, 2000)
#'
#' @importFrom wrapr stop_if_dot_args
#' @export
mean.D4Source <-
    function(x, chr, start, end, ..., track = NA) {
        stop_if_dot_args(substitute(list(...)), "mean()")
        x$mean(chr, start, end, track)
    }

#' For a set of genomic ranges in a `GRanges` object, adds the mean value over
#' the region as a metadata attribute for each genomic range, returning a new
#' `GRanges` object with updated metadata.
#'
#' @param d4source `D4Source` object
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @return New version of `granges` with region mean added as the `'mean'`
#' metadata attribute for each genomic range. If a track is specified, the
#' track name will be included in the metadata column name.
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' granges <- GRanges(seqnames = c('chr1'), ranges = IRanges(1000:2000), 
#'                    strand = c('+'))
#' granges_with_mean <- update_mean(d4source, granges)
#'
#' @importFrom GenomicRanges mcols seqnames start end ranges strand
#' @importFrom wrapr stop_if_dot_args
#' @export
update_mean <-
    function(d4source, granges, ..., track = NA) {
        stop_if_dot_args(substitute(list(...)), "update_mean()")
        fn <- function(chr, start, end, track) {
            mean.D4Source(
                d4source,
                chr = chr,
                start = start,
                end = end,
                track = track
            )
        }
        mcol_name <- if (is.na(track))
            "mean"
        else
            paste("mean", track, sep = "_")
        GenomicRanges::mcols(granges)[[mcol_name]] <-
            mapply(
                fn,
                chr = as.vector(seqnames(granges)),
                start = as.vector(start(ranges(granges))) - 1,
                end = as.vector(end(ranges(granges))),
                track = track,
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )

        granges

    }

#' Returns the median data value over the given region
#'
#' @param x `D4Source` object
#' @param na.rm Ignored; included in signature only to enable preservation of
#' `stats::median` in the namespace while implementing median for D4Source
#' @param ... Not used for values, forces later arguments to bind by name
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Median data value
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' median(d4source, chr = 'chr1', start = 1000, end = 2000)
#'
#' @importFrom wrapr stop_if_dot_args
#' @exportS3Method stats::median
median.D4Source <-
    function(x, na.rm = FALSE, ..., chr, start, end, track = NA) {
        stop_if_dot_args(substitute(list(...)), "median()")
        x$median(chr, start, end, track)
    }

#' For a set of genomic ranges in a `GRanges` object, adds the median value
#' over the region as a metadata attribute for each genomic range, returning a
#' new `GRanges` object with updated metadata.
#'
#' @param d4source `D4Source` object
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @return New version of `granges` with region median added as the `'median'`
#' metadata attribute for each genomic range. If a track is specified, the
#' track name will be included in the metadata column name.
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' granges <- GRanges(seqnames = c('chr1'), ranges = IRanges(1000:2000), 
#'                    strand = c('+'))
#' granges_with_median <- update_median(d4source, granges)
#'
#' @importFrom GenomicRanges mcols seqnames start end ranges strand
#' @importFrom wrapr stop_if_dot_args
#' @export
update_median <-
    function(d4source, granges, ..., track = NA) {
        stop_if_dot_args(substitute(list(...)), "update_median()")
        fn <- function(chr, start, end, track) {
            median.D4Source(
                d4source,
                chr = chr,
                start = start,
                end = end,
                track = track
            )
        }
        mcol_name <- if (is.na(track))
            "median"
        else
            paste("median", track, sep = "_")
        GenomicRanges::mcols(granges)[[mcol_name]] <-
            mapply(
                fn,
                chr = as.vector(seqnames(granges)),
                start = as.vector(start(ranges(granges))) - 1,
                end = as.vector(end(ranges(granges))),
                track = track,
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )

        granges

    }

#' Returns the data value at the given percentile for the given region
#'
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param percentile Percentile between 0 and 100. Example: percentile 0 is
#' the minimum data value, 50 is the median and 100 is the maximum.
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Data value at the given percentile
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' percentile(d4source, 'chr1', 1000, 2000, 50)
#'
#' @importFrom wrapr stop_if_dot_args
#' @export
percentile <-
    function(d4source, chr, start, end, percentile, ..., track = NA) {
        stop_if_dot_args(substitute(list(...)), "percentile()")
        d4source$percentile(chr, start, end, track, percentile)
    }

#' For a set of genomic ranges in a `GRanges` object, adds the percentile value
#' over the region as a metadata attribute for each genomic range, returning a
#' new `GRanges` object with updated metadata.
#'
#' @param d4source `D4Source` object
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param percentile Percentile between 0 and 100. Example: percentile 0 is the
#' minimum data value, 50 is the median and
#' 100 is the maximum.
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @return New version of `granges` with region percentile added as the
#' `'percentile_<percentile>'` metadata attribute for each genomic range. The
#' percentile requested (between 0 and 100) is included in the metadata column
#' name. If a track is specified, the track name will also be included in the
#' column name.
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' granges <- GRanges(seqnames = c('chr1'), ranges = IRanges(1000:2000), 
#'                    strand = c('+'))
#' granges_with_percentile <- update_percentile(d4source, granges, 50)
#'
#' @importFrom GenomicRanges mcols seqnames start end ranges strand
#' @importFrom wrapr stop_if_dot_args
#' @export
update_percentile <-
    function(d4source, granges, percentile, ..., track = NA) {
        stop_if_dot_args(substitute(list(...)), "update_percentile()")
        fn <- function(chr, start, end, track) {
            percentile(
                d4source,
                chr = chr,
                start = start,
                end = end,
                percentile = percentile,
                track = track
            )
        }
        mcol_name <- {
            if (is.na(track))
                paste("percentile", percentile, sep = "_")
            else
                paste("percentile", percentile, track, sep = "_")
        }
        GenomicRanges::mcols(granges)[[mcol_name]] <-
            mapply(
                fn,
                chr = as.vector(seqnames(granges)),
                start = as.vector(start(ranges(granges))) - 1,
                end = as.vector(end(ranges(granges))),
                track = track,
                SIMPLIFY = FALSE,
                USE.NAMES = FALSE
            )

        granges

    }


#' Bin across the given region, populating each bin with a summary of the
#' original data values in the bin region
#'
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param method Summary statistic for the data in each bin; one of 'mean' or
#' 'median'
#' @param ... Not used for values, forces later arguments to bind by name
#' @param bin_size Bin size
#' @param allow_bin_size_adjustment Only used for remote files: optimize bin
#' size for reading over a network
#' @param track Track name. If omitted or NA, queries the first track.
#' @return New query result of resampled data as a list
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' resample(d4source, 'chr1', 0, 50000000, 'mean')
#'
#' @export
resample <-
    function(
        d4source, chr, start, end, method, ..., bin_size = 1000,
        allow_bin_size_adjustment = TRUE, track = NA
    ) {
        stop_if_dot_args(substitute(list(...)), "resample()")
        query_result(
            d4source$resample(
                chr,
                start,
                end,
                track,
                method,
                bin_size,
                allow_bin_size_adjustment
            )
        )
    }
