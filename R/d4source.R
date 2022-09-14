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
      start = query_result_env$query()$left() + 1,
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

#' Prints information about a `D4Source` object
#'
#' @param d4source `D4Source` object
#' @return (Invisibly) the string that was printed
#' @import GenomeInfoDb
print.D4Source <- function(d4source) {
  seqs <- seqinfo(d4source)
  description <- paste(
    "<D4Source object>",
    paste("Source file:", get_source(d4source)),
    paste("Number of chromosomes:", length(seqs)),
    paste("Total chromosome length:", sum(seqlengths(seqs))),
    paste("Number of tracks:", length(get_tracks(d4source))),
    sep = "\n"
  )
  cat(description)
  invisible(description)
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

# Register a seqinfo() getter for `D4Source` class
setMethod("seqinfo", "D4Source", function(x) {
  chrs <- x$get_chroms()
  names <- vapply(chrs, function(c) c$name, "x")
  lengths <- vapply(chrs, function(c) c$size, 0)
  GenomeInfoDb::Seqinfo(names, seqlengths = lengths)
})

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
#' @param start 1-based start position
#' @param end 1-based inclusive end position
#' @param ... Not used for values, forces later arguments to bind by name
#' @param minus_strand If TRUE, vector of position scores will read in the
#' minus strand direction from `end` to `start`.
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
      d4source$query(chr, start - 1, end, track), minus_strand = minus_strand)
  }

#' Returns a `GRanges` object representing the data in a D4 track.
#' 
#' The returned `GRanges` contains ranges covering the entire requested region.
#' Runs of equal data values are collapsed to single ranges. Data values
#' are stored in a metadata column.
#' 
#' Note: if data values are continuous or vary a lot, for large regions this 
#' function may be slow and may return a very large `GRanges` object. In that 
#' case, consider working with smaller regions. Use this function when a 
#' `GRanges` object is required for other Bioconductor workflows.
#' 
#' @param d4source `D4Source` object
#' @param ... Not used for values, forces later arguments to bind by name
#' @param chr Chromosome/contig name. If omitted, queries the entire genome.
#' @param start 1-based start position. If omitted, queries the entire 
#' chromosome.
#' @param end 1-based inclusive end position. If omitted, queries the entire 
#' chromosome.
#' @param track Track name. If omitted or NA, queries the first track.
#' @param colname Name of metadata column in which to store 
#' @return `GRanges` object with metadata column storing data values 
#' @examples 
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' 
#' # Specific region
#' to_granges(d4source, chr = "chr1", start = 1000000, end = 1000000)
#' 
#' # Entire chromosome
#' to_granges(d4source, chr = "chr3")
#' 
#' # Entire genome
#' to_granges(d4source)
#' @importFrom wrapr stop_if_dot_args
#' @export
to_granges <- function(d4source, ..., chr = NA, start = NA, end = NA, track = NA, colname = "score") {
  
  if (is.na(chr)) {
    # No chromosome specified; query whole genome
    if (!is.na(start) || !is.na(end)) {
      warning("No chromosome provided. Ignoring start and/or end and querying whole genome.")
    }
    # Recursively apply to each chromosome and concatenate
    Reduce(
      append, 
      Map(function(c) to_granges(d4source, chr = c), seqnames(seqinfo(d4source))), 
      GRanges()
    )
  } else {
    
    # Check chromosome name is valid
    if (!(chr %in% seqnames(seqinfo(d4source)))) {
      stop(paste(chr, "not found"))
    }
    
    if (is.na(start) && is.na(end)) {
      # No coordinates provided; query entire chromosome
      start <- 1
      end <- seqlengths(seqinfo(d4source)[chr])[1]
    } 
    
    if (is.na(start) || is.na(end)) {
      # Insufficient coordinates provided; reset start and end to query entire chromosome
      warning("One of start/end is missing. Ignoring both and querying whole chromosome.")
      start <- 1
      end <- seqlengths(seqinfo(d4source)[chr])[1]
    }
    
    # Vector of data value per position
    values <- query(d4source, chr = chr, start = start, end = end, track = track)$results
    # Collapse runs of equal values
    r <- rle(values)
    r_lengths <- r$lengths
    r_values <- r$values
    # Build matrix of intervals with values
    mat <- matrix(nrow = length(r_lengths), ncol = 3)
    curr_start <- start
    for (curr_interval_idx in 1:nrow(mat)) {
      curr_end <- curr_start + r_lengths[curr_interval_idx] - 1
      mat[curr_interval_idx, 1] <- curr_start
      mat[curr_interval_idx, 2] <- curr_end
      mat[curr_interval_idx, 3] <- r_values[curr_interval_idx]
      curr_start <- curr_end + 1
    }
    # Convert to dataframe needed by makeGRangesFromDataFrame
    df <- data.frame(mat)
    colnames(df) <- c("start", "end", colname)
    df$chr <- chr
    # Build and return GRanges object
    GenomicRanges::makeGRangesFromDataFrame(
      df, keep.extra.columns = TRUE, ignore.strand = TRUE, seqinfo = seqinfo(d4source))
  }
}

#' For a set of genomic ranges in a `GRanges` object, extracts each region from
#' the D4 file and adds the query result as a metadata attribute for each
#' genomic range, returning a new `GRanges` object with updated metadata.
#'
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param d4source `D4Source` object
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @return New version of `granges` with query results added as the
#' `'query_results'` metadata attribute for each genomic range. For ranges on
#' the minus strand, the vector of position scores reads in the minus strand
#' direction. If a track is specified, the track name will be included in the
#' metadata column name.
#' @import GenomicRanges
#' @import IRanges
#' @importFrom wrapr stop_if_dot_args
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' granges <- GenomicRanges::GRanges(
#'              seqnames = c('chr1'), 
#'              ranges = IRanges::IRanges(1000:2000), 
#'              strand = c('+'))
#' granges_with_query <- update_with_query_results(granges, d4source)
#'
#' @export
update_with_query_results <-
  function(granges, d4source, ..., track = NA) {
    stop_if_dot_args(substitute(list(...)), "update_with_query_results()")
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
        start = as.vector(start(ranges(granges))),
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
#' @param start 1-based start position
#' @param end 1-based inclusive end position
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
    x$mean(chr, start - 1, end, track)
  }

#' For a set of genomic ranges in a `GRanges` object, adds the mean value over
#' the region as a metadata attribute for each genomic range, returning a new
#' `GRanges` object with updated metadata.
#'
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param d4source `D4Source` object
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @param colname Name of new metadata column to add. If omitted: the name of 
#' the new metadata column will be `"mean"` if no track was provided, or 
#' `"mean_<track>"` if a track was provided.
#' @return New version of `granges` with region mean added as a metadata 
#' attribute for each genomic range. 
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' granges <- GenomicRanges::GRanges(
#'              seqnames = c('chr1'), 
#'              ranges = IRanges::IRanges(1000:2000), 
#'              strand = c('+'))
#' granges_with_mean <- update_with_mean(granges, d4source)
#'
#' @importFrom GenomicRanges mcols seqnames start end ranges strand
#' @importFrom wrapr stop_if_dot_args
#' @export
update_with_mean <-
  function(granges, d4source, ..., track = NA, colname = NA) {
    stop_if_dot_args(substitute(list(...)), "update_with_mean()")
    fn <- function(chr, start, end, track) {
      mean.D4Source(
        d4source,
        chr = chr,
        start = start,
        end = end,
        track = track
      )
    }
    mcol_name <- if (!is.na(colname)) colname 
    else if (is.na(track))
      "mean"
    else
      paste("mean", track, sep = "_")
    GenomicRanges::mcols(granges)[[mcol_name]] <-
      mapply(
        fn,
        chr = as.vector(seqnames(granges)),
        start = as.vector(start(ranges(granges))),
        end = as.vector(end(ranges(granges))),
        track = track,
        SIMPLIFY = TRUE,
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
#' @param start 1-based start position
#' @param end 1-based inclusive end position
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
    x$median(chr, start - 1, end, track)
  }

#' For a set of genomic ranges in a `GRanges` object, adds the median value
#' over the region as a metadata attribute for each genomic range, returning a
#' new `GRanges` object with updated metadata.
#'
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param d4source `D4Source` object
#' @param ... Not used for values, forces later arguments to bind by name
#' @param track Track name. If omitted or NA, queries the first track.
#' @param colname Name of new metadata column to add. If omitted: the name of 
#' the new metadata column will be `"median"` if no track was provided, or 
#' `"median_<track>"` if a track was provided.
#' @return New version of `granges` with region median added as a metadata 
#' attribute for each genomic range. 
#' @examples
#' d4source <- D4Source(system.file('extdata', 'example.d4', package = 'rd4'))
#' granges <- GenomicRanges::GRanges(
#'              seqnames = c('chr1'), 
#'              ranges = IRanges::IRanges(1000:2000), 
#'              strand = c('+'))
#' granges_with_median <- update_with_median(granges, d4source)
#'
#' @importFrom GenomicRanges mcols seqnames start end ranges strand
#' @importFrom wrapr stop_if_dot_args
#' @export
update_with_median <-
  function(granges, d4source, ..., track = NA, colname = NA) {
    stop_if_dot_args(substitute(list(...)), "update_with_median()")
    fn <- function(chr, start, end, track) {
      median.D4Source(
        d4source,
        chr = chr,
        start = start,
        end = end,
        track = track
      )
    }
    mcol_name <- if(!is.na(colname)) colname
    else if (is.na(track))
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
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
    
    granges
    
  }

#' Returns the data value at the given percentile for the given region
#'
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 1-based start position
#' @param end 1-based inclusive end position
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
    d4source$percentile(chr, start - 1, end, track, percentile)
  }

#' For a set of genomic ranges in a `GRanges` object, adds the percentile value
#' over the region as a metadata attribute for each genomic range, returning a
#' new `GRanges` object with updated metadata.
#'
#' @param granges `GRanges` object containing genomic ranges of interest
#' @param d4source `D4Source` object
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
#' granges <- GenomicRanges::GRanges(
#'              seqnames = c('chr1'), 
#'              ranges = IRanges::IRanges(1000:2000), 
#'              strand = c('+'))
#' granges_with_percentile <- update_with_percentile(granges, d4source, 50)
#'
#' @importFrom GenomicRanges mcols seqnames start end ranges strand
#' @importFrom wrapr stop_if_dot_args
#' @export
update_with_percentile <-
  function(granges, d4source, percentile, ..., track = NA) {
    stop_if_dot_args(substitute(list(...)), "update_with_percentile()")
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
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
    
    granges
    
  }


#' Bin across the given region, populating each bin with a summary of the
#' original data values in the bin region
#'
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 1-based start position
#' @param end 1-based inclusive end position
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
        start - 1,
        end,
        track,
        method,
        bin_size,
        allow_bin_size_adjustment
      )
    )
  }

