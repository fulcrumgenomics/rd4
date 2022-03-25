//! R bindings for the [D4](https://github.com/38/d4-format) file format.
//!
//! See R documentation for how to use this.

mod d4_reader;
use d4::Chrom;
use d4_reader::{http_d4_reader::HttpD4Reader, local_d4_reader::LocalD4Reader, D4Reader};
use extendr_api::prelude::*;
use serde::{Deserialize, Serialize};

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
    /// @export
    fn query(&self, chr: String, left: u32, right: u32, track: Option<String>) -> QueryResult {
        self.inner.query(chr, left, right, track)
    }
}

/// Helper struct to hold onto the result of a query, and the context of the query.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryResult {
    /// The depths from the query
    result: Vec<i32>,
    /// The query used to generate this result
    query: Query,
    /// The d4 source being queried
    source: String,
    /// The track of the d4 source, if specified
    track: Option<String>,
}

impl QueryResult {
    /// Create a new [`QueryResult`]
    fn new(query: Query, d4_file_path: String, d4_track: Option<String>, result: Vec<i32>) -> Self {
        Self { result, query, source: d4_file_path, track: d4_track }
    }
}

/// The object returned after querying a [`D4Source`]
/// @export
#[extendr]
impl QueryResult {
    /// Get the results of the query
    ///
    /// @export
    fn results(&self) -> &[i32] {
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

/// Wrapper type for converting a [`Chrom`] into an R List.
struct ChromWrapper(Chrom);
impl From<ChromWrapper> for Robj {
    fn from(c: ChromWrapper) -> Self {
        // The cast o i32 forces this to be a "Real" integer in R, otherwise it is treated as a real (aka float)
        r!(list!(name = c.0.name, size = c.0.size as i32))
    }
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rd4;
    impl D4Source;
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
            let tracks = file.get_tracks();
            eprintln!("Length is {:?}", tracks.len());
            assert_eq!(file.get_tracks().elt(0), "");
            assert_eq!(file.get_chroms().len(), 1);
            assert_eq!(file.get_chroms().elt(0)?.dollar("name")?, r!("chr1"));
            assert_eq!(file.get_chroms().elt(0)?.dollar("size")?, r!(1000_i32));
            let r = file.query(String::from("chr1"), 12, 22, None);
            assert_eq!(r.results(), &[0, 0, 0, 200, 0, 0, 0, 0, 100, 0]);
            assert_eq!(r.source(), data.to_string_lossy());
            assert_eq!(r.track(), String::from(""));
            let q = r.query();
            assert_eq!(q, Query::new(String::from("chr1"), 12, 22));
        }
    }
}
