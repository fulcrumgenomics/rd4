# Convert a Rust QueryResultEnv to an R list
query_result <- function(query_result_env) {
  list(
    results = query_result_env$results(),
    query = list(
      chr = query_result_env$query()$chr(),
      start = query_result_env$query()$left(),
      end = query_result_env$query()$right()
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
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' methods(class="D4Source")
#' 
#' @export
D4Source <- function(path) {
  d4 <- D4SourceEnv$new(path)
  class(d4) <- c("D4Source", class(d4))
  d4
}

#' @export
get_source <- function(obj, ...) {UseMethod("get_source")}

#' Returns the path to the source D4 file
#' 
#' @return The path to the source D4 file
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' get_source(d4source)
#' 
#' @export
get_source.D4Source <- function(d4source) {d4source$get_source()}

#' @export
get_chroms <- function(obj, ...) {UseMethod("get_chroms")}

#' Returns the list of chromosomes in the source D4 file
#' 
#' @return Nested list; each element represents a chromosome and stores the chromosome name and length
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' get_chroms(d4source)
#' 
#' @export
get_chroms.D4Source <- function(d4source) {d4source$get_chroms()}

#' @export
get_tracks <- function(obj, ...) {UseMethod("get_tracks")}

#' Returns the names of tracks in the source D4 file
#' 
#' @return Vector of track names
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' get_tracks(d4source)
#' 
#' @export
get_tracks.D4Source <- function(d4source) {d4source$get_tracks()}

#' @export
query <- function(obj, ...) {UseMethod("query")}

#' Extracts the given region from the D4 file
#' 
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Query result as a list
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' query_result <- query(d4source, "chr1", 1000, 2000)
#' 
#' @export
query.D4Source <- function(d4source, chr, start, end, track = NA) {
  query_result(d4source$query(chr, start, end, track))
}

#' Returns the mean data value over the given region
#' 
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Mean data value
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' mean(d4source, "chr1", 1000, 2000)
#' 
#' @export
mean.D4Source <- function(d4source, chr, start, end, track = NA) {d4source$mean(chr, start, end, track)}

#' Returns the median data value over the given region
#' 
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Mean data value
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' median(d4source, "chr1", 1000, 2000)
#' 
#' @export
median.D4Source <- function(d4source, chr, start, end, track = NA) {d4source$median(chr, start, end, track)}

#' @export
percentile <- function(obj, ...) {UseMethod("percentile")}

#' Returns the data value at the given percentile for the given region
#' 
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param percentile Percentile between 0 and 100. Example: percentile 0 is the minimum data value, 50 is the median and 
#' 100 is the maximum.
#' @param track Track name. If omitted or NA, queries the first track.
#' @return Data value at the given percentile
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' percentile(d4source, "chr1", 1000, 2000, 50)
#' 
#' @export
percentile.D4Source <- function(d4source, chr, start, end, percentile, track = NA) {
  d4source$percentile(chr, start, end, track, percentile)
}

#' @export
resample <- function(obj, ...) {UseMethod("resample")}

#' Bin across the given region, populating each bin with a summary of the original data values in the bin region
#' 
#' @param d4source `D4Source` object
#' @param chr Chromosome/contig name
#' @param start 0-based start position
#' @param end 0-based exclusive end position
#' @param method Summary statistic for the data in each bin; one of 'mean' or 'median'
#' @param bin_size Bin size
#' @param allow_bin_size_adjustment Only used for remote files: optimize bin size for reading over a network
#' @param track Track name. If omitted or NA, queries the first track.
#' @return New query result of resampled data as a list
#' @examples
#' d4source <- D4Source(system.file("extdata", "example.d4", package="rd4"))
#' resample(d4source, "chr1", 0, 50000000, "mean")
#' 
#' @export
resample.D4Source <- function(
  d4source, 
  chr, 
  start, 
  end, 
  method, 
  bin_size = 1000, 
  allow_bin_size_adjustment = TRUE, 
  track = NA
) {
  query_result(d4source$resample(chr, start, end, track, method, bin_size, allow_bin_size_adjustment))
}

