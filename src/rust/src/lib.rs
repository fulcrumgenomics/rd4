//! R bindings for the [D4](https://github.com/38/d4-format) file format.
//!
//! See R documentation for how to use this.

// TODO: broader "range" parsing more like python api where you can give just a chr and it will get the length for you

mod d4_reader;
use d4::Chrom;
use d4_reader::{http_d4_reader::HttpD4Reader, local_d4_reader::LocalD4Reader, D4Reader};
use extendr_api::prelude::*;
use serde::{Deserialize, Serialize};

/// Helper fuction to convert an [`Robj`] to an i32.
fn robj_to_i32(obj: &Robj) -> Option<i32> {
    if obj.is_na() {
        return None;
    }
    // There does not appear to be a transparent conversion from RObj to i32 that just works from the R side with no annoatation.
    // Since R treats numbers as f64 by default, that's what we accept here and then cast to i32.
    if let Ok(value) = <Option<f64>>::from_robj(obj) {
        return value.map(|v| v as i32);
    }
    if let Ok(value) = <Option<i32>>::from_robj(obj) {
        return value;
    }
    panic!("Expected numeric object, found: {:?}", obj);
}

/// The R object for reading from a D4 source.
///
/// **Note**: the source path usage in these bindings differs from both d4 and pyd4 in that this library
/// expects the "track" to be passed in separately when performing a query, instead of creating a "D4File"
/// that implicitly opens a specific track that is appended to the path.
#[derive(Debug)]
struct D4Source {
    source: String,
    inner: Box<dyn D4Reader>,
}

/// D4Source class
///
/// This is the primary way of interacting with D4 sources.
///
/// @examples
/// d4 <- D4File$new(path)
/// chroms <- d4$get_chroms()
/// @export
#[extendr(use_try_from = true)]
impl D4Source {
    /// Method for making a new [`D4Source`]  object.
    ///
    /// **Note**: This expects a plain path, with no track appended to the end.
    /// @export
    fn new(source: &str) -> Self {
        let inner: Box<dyn D4Reader> =
            if source.starts_with("http://") || source.starts_with("https://") {
                Box::new(HttpD4Reader::new(source.to_string()))
            } else {
                Box::new(LocalD4Reader::new(source.to_string()))
            };
        Self { source: source.to_string(), inner }
    }

    /// Method for getting the original source path used to create this object
    ///
    /// @export
    fn get_source(&self) -> String {
        self.source.clone()
    }

    /// List all chromosomes in the [`D4Source`]
    ///
    /// @export
    fn get_chroms(&self) -> List {
        List::from_values(self.inner.get_chroms().into_iter().map(ChromWrapper))
    }

    /// @export
    fn custom(&self) -> String {
        String::from("CUSTOM")
    }

    /// List all tracks in the [`D4Source`]
    /// @export
    fn get_tracks(&self) -> Strings {
        Strings::from_values(self.inner.get_tracks())
    }

    /// Query the [`D4Source`] for the depths over the specified region
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if NA, the first track in the source will be queried
    ///
    /// @export
    fn query(&self, chr: String, left: u32, right: u32, track: Option<String>) -> QueryResult {
        self.inner.query(chr, left, right, track)
    }

    /// Compute the mean depth of the given region.
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if NA, the first track in the source will be queried
    ///
    /// @export
    fn mean(&self, chr: String, left: u32, right: u32, track: Option<String>) -> f64 {
        let query = Query::new(chr, left, right);
        self.inner.mean(&[query], track.as_deref())[0]
    }

    /// Compute the median depth of the given region.
    ///
    /// **Note** this will fail on remote files and is only supported for local files.
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if NA, the first track in the source will be queried
    ///
    /// @export
    fn median(&self, chr: String, left: u32, right: u32, track: Option<String>) -> f64 {
        let query = Query::new(chr, left, right);
        if let Some(reader) = self.inner.as_any().downcast_ref::<LocalD4Reader>() {
            reader.percentile(&[query], track.as_deref(), 50.0)[0]
        } else {
            panic!("Median operation is not supported on remote D4 sources.")
        }
    }

    /// Returns the value histogram for the given region.
    ///
    /// **Note** this operation is only supported on local D4 sources.
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if NA, the first track in the source will be queried
    /// - `min` - the smallest bucket of the histogram
    /// - `max` - the largest bucket of the histogram, 1024 if None
    ///
    /// @export
    fn histogram(
        &self,
        chr: String,
        left: u32,
        right: u32,
        track: Option<String>,
        min: i32,
        max: Robj,
    ) -> Histogram {
        let query = Query::new(chr, left, right);
        let max = robj_to_i32(&max);
        if let Some(reader) = self.inner.as_any().downcast_ref::<LocalD4Reader>() {
            reader
                .histogram(&[query], min, max.unwrap_or(1024), track.as_deref())
                .into_iter()
                .next()
                .unwrap()
        } else {
            panic!("Histogram operation is not supported on remote D4 sources.")
        }
    }

    /// Returns the percentile value in the given region.
    ///
    /// **Note** this operation is only supported on local D4 sources.
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if NA, the first track in the source will be queried
    /// - `percentile` - the percentile value
    ///
    /// @export
    fn percentile(
        &self,
        chr: String,
        left: u32,
        right: u32,
        track: Option<String>,
        percentile: f64,
    ) -> f64 {
        let query = Query::new(chr, left, right);
        if let Some(reader) = self.inner.as_any().downcast_ref::<LocalD4Reader>() {
            reader.percentile(&[query], track.as_deref(), percentile)[0]
        } else {
            panic!("Histogram operation is not supported on remote D4 sources.")
        }
    }

    /// Re-sample the given region.
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if NA, the first track in the source will be queried
    /// - `method` - the method to use for resampling, mean and median are supported. Only mean is supported for remote files
    /// - `bin_size` - optional bin_size argument sets the width of bins to look at for re-sampling
    /// - `allow_bin_size_adjustment` - only used for remote files, this will optimize the bin size for reading over a network
    ///
    /// @export
    #[allow(clippy::too_many_arguments)]
    fn resample(
        &self,
        chr: String,
        left: u32,
        right: u32,
        track: Option<String>,
        method: String,
        bin_size: Robj,
        allow_bin_size_adjustment: Option<bool>,
    ) -> QueryResult {
        let bin_size = robj_to_i32(&bin_size);
        let bin_size = self.inner.adjust_bin_size(
            bin_size.unwrap_or(1_000) as u32,
            allow_bin_size_adjustment.unwrap_or(true),
        );
        let query = Query::new(chr, left, right);
        assert!(bin_size > 0, "Invalid bin_size argument. Must be greater than 0.");

        // Split the given region into bin size chunks
        let mut splitted = vec![];
        let mut begin = query.left;
        let end = query.right;
        while begin < end {
            let bin_end = std::cmp::min(begin + bin_size, end);
            splitted.push(Query::new(query.chr.clone(), begin, bin_end));
            begin = bin_end;
        }

        let values = if method == "mean" {
            self.inner.mean(&splitted, track.as_deref())
        } else if method == "median" {
            if let Some(reader) = self.inner.as_any().downcast_ref::<LocalD4Reader>() {
                reader.percentile(&splitted, track.as_deref(), 50.0)
            } else {
                panic!("Re-sampling by median is only supported with local D4 Files")
            }
        } else {
            panic!(
                "{}, is an unsupported resampling method. Please choose one of [mean, median]",
                method
            );
        };

        QueryResult::new(query, self.get_source(), track, values, Some(bin_size as i32))
    }
}

/// Helper struct to hold onto the result of a query, and the context of the query.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryResult {
    /// The depths from the query
    result: Vec<f64>,
    /// The query used to generate this result
    query: Query,
    /// The d4 source being queried
    source: String,
    /// The track of the d4 source, if specified
    track: Option<String>,
    /// The bin size used for this query (if resampled)
    bin_size: Option<i32>,
}

impl QueryResult {
    /// Create a new [`QueryResult`]
    fn new(
        query: Query,
        d4_file_path: String,
        d4_track: Option<String>,
        result: Vec<f64>,
        bin_size: Option<i32>,
    ) -> Self {
        Self { result, query, source: d4_file_path, track: d4_track, bin_size }
    }
}

/// The object returned after querying a [`D4Source`]
/// @export
#[extendr]
impl QueryResult {
    /// Get the results of the query
    ///
    /// @export
    fn results(&self) -> &[f64] {
        &self.result
    }

    /// Return the [`Query`] used to get this result
    ///
    /// @export
    fn query(&self) -> Query {
        self.query.clone()
    }

    /// Get the D4 file that was queried
    ///
    /// @export
    fn source(&self) -> String {
        self.source.clone()
    }

    /// Get the D4 track that was queried
    ///
    /// Returns and empty string if none was specified
    ///
    /// @export
    fn track(&self) -> String {
        self.track.clone().unwrap_or_else(|| String::from(""))
    }

    /// Get the bin size used for this query.
    ///
    /// This is only useful if the query originated from `resample`, otherwise it will be NA.
    ///
    /// @export
    fn bin_size(&self) -> Option<i32> {
        self.bin_size
    }
}

/// A [`Query`] holds onto user input `chr`, `left`, and `right` specifiers.
#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Serialize, Deserialize)]
struct Query {
    /// The specified chromosome
    chr: String,
    /// The inclusive start of the query
    left: u32,
    /// The exclusive stop of the query
    right: u32,
}

impl Query {
    /// Create a new query
    fn new(chr: String, left: u32, right: u32) -> Self {
        Self { chr, left, right }
    }
}

/// A context object that holds onto the original query parameters
/// @export
#[extendr]
impl Query {
    /// Get the chromosome specified by this query
    ///
    /// @export
    fn chr(&self) -> String {
        self.chr.clone()
    }

    /// Get the inclusive "left" position specified by this query
    ///
    /// @export
    fn left(&self) -> u32 {
        self.left
    }

    /// Get the inclusive "right" position specified by this query
    ///
    /// @export
    fn right(&self) -> u32 {
        self.right
    }
}

impl<'a> From<&'a Query> for (&'a str, u32, u32) {
    fn from(val: &'a Query) -> Self {
        (val.chr.as_ref(), val.left, val.right)
    }
}

/// Wrapper type for converting a [`Chrom`] into an R List.
struct ChromWrapper(Chrom);
impl From<ChromWrapper> for Robj {
    fn from(c: ChromWrapper) -> Self {
        // The cast o i32 forces this to be a "Real" integer in R, otherwise it is treated as a real (aka float)
        r!(list!(name = c.0.name, size = c.0.size as i32))
    }
}

/// Represents a Histogram of counts of the values between the `min` and `max` used when creating the histogram
///
///
#[derive(Debug, Clone)]
struct Histogram {
    below: u32,
    above: u32,
    first_value: i32,
    prefix_sum: Vec<u32>,
}

impl Histogram {
    /// Create a new histogram object with a vec of values and counts.
    fn new(mut values: Vec<(i32, u32)>, below: u32, above: u32) -> Self {
        values.sort_unstable(); // TODO: I don't think this sort is needed, the "values" is the (value, count), and should be in value order I think?
        let first_value = values[0].0;
        let mut prefix_sum = vec![below];
        let mut current_value = first_value;
        for (value, count) in values.into_iter() {
            while current_value < value - 1 {
                current_value += 1;
                prefix_sum.push(*prefix_sum.last().unwrap());
            }
            current_value += 1;
            prefix_sum.push(prefix_sum.last().unwrap() + count)
        }

        Self { below, above, first_value, prefix_sum }
    }
}

/// The [`Histogram`] API.
/// @export
#[extendr]
impl Histogram {
    /// Get the counts of a value in the histogram.
    ///
    /// If the value is outside the bounds specified when creating the histogram, `0` will be returned.
    ///
    /// @export
    fn value_count(&self, value: i32) -> i32 {
        // TODO: should the "value" be multiplied by the denominator, and be f64?
        if value < self.first_value || self.first_value + (self.prefix_sum.len() as i32) - 1 < value
        {
            0
        } else {
            let idx = (value - self.first_value + 1) as usize;
            (self.prefix_sum[idx] - self.prefix_sum[idx - 1]) as i32
        }
    }

    /// Total number of data points, including points above and below the specified histogram range.
    ///
    /// @export
    fn total_count(&self) -> i32 {
        (self.prefix_sum.last().unwrap() + self.above) as i32
    }

    /// The fraction represented by the the number of counts of the specified value out of the total counts (including those above and below the specified histogram range).
    ///
    /// @export
    fn value_fraction(&self, value: i32) -> f64 {
        // TODO: same question about value and denominator
        self.value_count(value) as f64 / self.total_count() as f64
    }

    /// The fraction of values below the specified value (includes counts below the min of the specified range) out of the total counts (including those above and below the specified range).
    /// @export
    fn fraction_below(&self, value: i32) -> f64 {
        if value < self.first_value || self.first_value + self.prefix_sum.len() as i32 - 1 < value {
            0.0
        } else {
            let idx = (value - self.first_value + 1) as usize;
            self.prefix_sum[idx] as f64 / self.total_count() as f64
        }
    }

    /// Get the mean from this histogram.
    ///
    /// @export
    fn mean(&self) -> f64 {
        let mut current_value = self.first_value as u32;
        let mut current_sum = self.prefix_sum[0];
        let mut total = 0;
        for value in self.prefix_sum.iter().skip(1) {
            let current_count = value - current_sum;
            total += current_count * current_value;
            current_value += 1;
            current_sum = *value;
        }
        total as f64 / self.total_count() as f64
    }

    /// Get the n-th percentile value from the histogram.
    ///
    /// **Note**, if the percentile falls outside the range specified for the histogram
    /// `min` and `max` then the `min` or `max`, respectively, will be returned.
    ///
    /// @export
    fn percentile(&self, percentile: f64) -> i32 {
        let total_count = self.total_count() as f64;
        let mut value = self.first_value;
        for count in self.prefix_sum.iter().skip(1) {
            if *count as f64 * 100.0 / total_count > percentile {
                return value;
            }
            value += 1;
        }
        0
    }

    /// Compute the median value of the histogram (i.e. the median depth).
    ///
    /// @export
    fn median(&self) -> i32 {
        self.percentile(50.0)
    }

    /// Compute the standard deviation of this histogram.
    ///
    /// @export
    fn std(&self) -> f64 {
        let mut current_value = self.first_value as u32;
        let mut current_sum = self.prefix_sum[0];
        let mut sum = 0;
        let mut square_sum = 0;
        for value in self.prefix_sum.iter().skip(1) {
            let current_count = value - current_sum;
            sum += current_count * current_value;
            square_sum += current_count * current_value * current_value;
            current_value += 1;
            current_sum = *value
        }

        let ex = sum as f64 / self.total_count() as f64;
        let esx = square_sum as f64 / self.total_count() as f64;
        (esx - ex * ex).sqrt()
    }
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rd4;
    impl D4Source;
    impl Histogram;
    impl Query;
    impl QueryResult;
}

#[cfg(test)]
mod test {
    use std::path::{Path, PathBuf};

    use d4::{
        ptab::PTablePartitionWriter, stab::SecondaryTablePartWriter, Chrom, D4FileBuilder,
        D4FileWriter,
    };
    use d4tools::make_dictionary;
    use extendr_api::prelude::*;
    use tempfile::TempDir;

    use crate::{D4Source, Query};

    fn create_d4_file<P: AsRef<Path>>(dir: P) -> PathBuf {
        let values = vec![(10, 100), (15, 200), (20, 100), (50, 1000)];

        let filename = dir.as_ref().join("test.d4");
        let mut builder = D4FileBuilder::new(&filename);
        builder.set_dictionary(make_dictionary(None, None).unwrap());
        builder.append_chrom(std::iter::once(Chrom { name: String::from("chr1"), size: 1000 }));
        let mut writer: D4FileWriter = builder.create().unwrap();
        // If there is a better way to get access to the inner tables I haven't found it
        let partitions = writer.parallel_parts(Some(10_000)).unwrap();
        assert_eq!(partitions.len(), 1);

        let (mut primary_table, mut secondary_table) = partitions.into_iter().next().unwrap();
        let mut p_encoder = primary_table.make_encoder();
        for (pos, depth) in values {
            if !p_encoder.encode(pos, depth) {
                secondary_table.encode(pos as u32, depth).unwrap();
            }
        }
        secondary_table.flush().unwrap();
        secondary_table.finish().unwrap();

        filename
    }

    #[test]
    fn test_simple_case() {
        test! {
            let tempdir = TempDir::new().unwrap();
            let data = create_d4_file(tempdir.path());

            let file = D4Source::new(data.to_str().unwrap());

            assert_eq!(file.get_tracks().elt(0), "");
            assert_eq!(file.get_chroms().len(), 1);
            assert_eq!(file.get_chroms().elt(0)?.dollar("name")?, r!("chr1"));
            assert_eq!(file.get_chroms().elt(0)?.dollar("size")?, r!(1000_i32));

            let r = file.query(String::from("chr1"), 12, 22, None);
            assert_eq!(r.results(), &[0., 0., 0., 200., 0., 0., 0., 0., 100., 0.]);
            assert_eq!(r.source(), data.to_string_lossy());
            assert_eq!(r.track(), String::from(""));

            assert_eq!(file.mean(String::from("chr1"), 12, 22, None), 30.0);
            assert_eq!(file.median(String::from("chr1"), 12, 22, None), 0.0);
            assert_eq!(file.percentile(String::from("chr1"), 12, 22, None, 50.0), 0.0);

            let resampled = file.resample(String::from("chr1"), 12, 22, None, String::from("mean"), r!(i32::na()), None);
            assert_eq!(resampled.result.len(), 1);

            let resampled = file.resample(String::from("chr1"), 12, 22, None, String::from("median"), r!(i32::na()), None);
            assert_eq!(resampled.result.len(), 1);

            let resampled = file.resample(String::from("chr1"), 0, 1000, None, String::from("mean"), r!(10), None);
            assert_eq!(resampled.result.len(), 100);
            let resampled = file.resample(String::from("chr1"), 0, 1000, None, String::from("median"), r!(10), None);
            assert_eq!(resampled.result.len(), 100);
        }
    }

    #[test]
    fn test_histogram_partial_range() {
        test! {
            let tempdir = TempDir::new().unwrap();
            let data = create_d4_file(tempdir.path());

            let file = D4Source::new(data.to_str().unwrap());
            let hist = file.histogram(String::from("chr1"), 0, 1000, None, 99, r!(200));
            assert_eq!(hist.below, 996);
            assert_eq!(hist.above, 2);
            assert_eq!(hist.mean(), 0.2);
            assert_eq!(hist.median(), 99);
            assert_eq!(hist.percentile(99.79), 100); // The two counts at depth 100 push the total observed seen points to 998
            assert_eq!(hist.fraction_below(199), 0.998);
            assert_eq!(hist.value_count(100), 2);
            assert_eq!(hist.value_fraction(100), 0.002);
            assert_eq!(hist.std(), 4.467_661_580_737_736);
        }
    }

    #[test]
    fn test_histogram_full_range() {
        test! {
            let tempdir = TempDir::new().unwrap();
            let data = create_d4_file(tempdir.path());

            let file = D4Source::new(data.to_str().unwrap());
            let hist = file.histogram(String::from("chr1"), 0, 1000, None, 0, r!(i32::na()));
            assert_eq!(hist.below, 0);
            assert_eq!(hist.above, 0);
            assert_eq!(hist.mean(), 1.4);
            assert_eq!(hist.median(), 0);
            assert_eq!(hist.percentile(99.79), 100); // The two counts at depth 100 push the total observed seen points to 998
            assert_eq!(hist.fraction_below(199), 0.998);
            assert_eq!(hist.value_count(100), 2);
            assert_eq!(hist.value_fraction(100), 0.002);
            assert_eq!(hist.std(), 32.527_526_804_231_52);
        }
    }
}
