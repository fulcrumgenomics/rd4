use extendr_api::prelude::*;

use d4::ptab::DecodeResult;
use d4::stab::SecondaryTablePartReader;
use d4::Chrom;

enum ReaderWrapper {
    LocalReader(d4::D4TrackReader),
    RemoteReader(d4::ssio::D4TrackReader<d4::ssio::http::HttpReader>),
}

impl ReaderWrapper {
    fn open(path: &str, is_remote: bool, track: Option<&str>) -> ReaderWrapper {
        if is_remote {
            let conn = d4::ssio::http::HttpReader::new(path).unwrap();
            let reader = d4::ssio::D4TrackReader::from_reader(conn, track).expect("todo");
            Self::RemoteReader(reader)
        } else {
            // TODO: D4TrackReader::open splits the path on `:` and takes anything to the right as a track specifier
            // There should be a way to create a reader like remote with the optional track name
            let track_spec =
                if let Some(track) = track { [path, track].join(":") } else { path.to_string() };
            let local_reader = d4::D4TrackReader::open(&track_spec).expect("todo");
            Self::LocalReader(local_reader)
        }
    }

    fn get_chroms(&self) -> &[Chrom] {
        match self {
            Self::LocalReader(local) => local.header().chrom_list(),
            Self::RemoteReader(remote) => remote.chrom_list(),
        }
    }
}

/// The R object for reading a D4 file
#[derive(Debug)]
struct D4File {
    path: String,
    is_remote: bool,
    track: Option<String>,
}

impl D4File {
    pub(crate) fn open(&self) -> ReaderWrapper {
        ReaderWrapper::open(self.path.as_str(), self.is_remote, self.track.as_deref())
    }

    /// Parse the path to try to extract a track
    fn parse_path(path: &str) -> (String, Option<String>) {
        if path.starts_with("http://") || path.starts_with("https://") {
            if let Some(split_pos) = path.rfind('#') {
                (path[..split_pos].to_string(), Some(path[split_pos + 1..].to_string()))
            } else {
                (path.to_string(), None)
            }
        } else if let Some(split_pos) = path.rfind(':') {
            (path[..split_pos].to_string(), Some(path[split_pos + 1..].to_string()))
        } else {
            (path.to_string(), None)
        }
    }
}

/// D4File class
/// @examples
/// d4 <- D4File$new(path)
/// @export
#[extendr]
impl D4File {
    /// Method for making a new object
    /// @export
    // TODO: I'm sticking with the same API philosohy as pyd4, which allows specifying a track in the file path
    fn new(path: &str) -> Self {
        let is_remote = path.starts_with("http://") || path.starts_with("https://");
        let (path, track) = Self::parse_path(path);
        Self { path, is_remote, track }
    }

    /// Method for getting the path
    /// @export
    fn get_path(&self) -> String {
        self.path.clone()
    }

    /// List all chromosomes in the [`D4File`]
    /// @export
    fn chroms(&self) -> Vec<String> {
        self.open().get_chroms().iter().map(|x| x.name.clone()).collect()
    }

    /// List all tracks in the [`D4File`]
    /// @export
    fn list_tracks(&self) -> Vec<String> {
        let mut tracks = vec![];
        d4::find_tracks_in_file(&self.path, |_| true, &mut tracks)
            .expect("Failed to extract tracks from file.");
        tracks.into_iter().map(|x| x.to_string_lossy().to_string()).collect()
    }

    /// Get values in region
    ///
    /// To get the values for a specific track, the `D4File` must be opened with a track_spec on the path.
    /// i.e. `/path/to/file.d4:<track_name>`
    /// @export
    // TODO - this is only returning data for a single track
    // TODO - Add query and file info to return type
    // TODO - I've used the same naming as pyd4 with left and right
    fn query(&self, chr: &str, left: u32, right: u32) -> QueryResult {
        let result = match self.open() {
            ReaderWrapper::RemoteReader(mut inner) => {
                let iter: Box<dyn Iterator<Item = i32>> =
                    Box::new(inner.get_view(chr, left, right).expect("todo").map(|res| {
                        if let Ok((_, value)) = res {
                            value
                        } else {
                            0
                        }
                    }));
                iter
            }
            ReaderWrapper::LocalReader(mut inner) => {
                let partition = inner.split(None).expect("todo");
                let chr = chr.to_string();
                let iter: Box<dyn Iterator<Item = i32>> =
                    Box::new(partition.into_iter().flat_map(move |(mut ptab, mut stab)| {
                        let (part_chr, begin, end) = ptab.region();
                        let part_chr = part_chr.to_string();
                        let pd = ptab.to_codec();
                        (if part_chr == chr { left.max(begin)..right.min(end) } else { 0..0 }).map(
                            move |pos| match pd.decode(pos as usize) {
                                DecodeResult::Definitely(value) => value,
                                DecodeResult::Maybe(value) => {
                                    if let Some(st_value) = stab.decode(pos) {
                                        st_value
                                    } else {
                                        value
                                    }
                                }
                            },
                        )
                    }));
                iter
            }
        }
        .collect();
        QueryResult::new(
            Query::new(chr.to_string(), left, right),
            self.path.clone(),
            self.track.clone(),
            result,
        )
    }
}

/// Helper struct to hold onto the result of a query, and the context of the query.
struct QueryResult {
    /// The depths from the query
    result: Vec<i32>,
    /// The query used to generate this result
    query: Query,
    /// The d4 file being queried
    d4_file_path: String,
    /// The track of the d4 file, if specified
    d4_track: Option<String>,
}

impl QueryResult {
    // TODO: should the result type be more generic?
    fn new(query: Query, d4_file_path: String, d4_track: Option<String>, result: Vec<i32>) -> Self {
        Self { result, query, d4_file_path, d4_track }
    }
}

/// The object returned after querying a [`D4File`]
/// @export
#[extendr]
impl QueryResult {
    /// Get the results of the query
    /// @export
    pub fn results(&self) -> &[i32] {
        &self.result
    }

    /// Return the query used to get this result: (chr, left, right)
    /// @export
    // TODO - return a tuple? make a query obj?
    pub fn query(&self) -> &Query {
        &self.query
    }

    /// Get the D4 file that was queried
    /// @export
    pub fn d4_file(&self) -> &str {
        &self.d4_file_path
    }

    /// Get the D4 track that was queried
    ///
    /// Returns and empty string if none was specified
    /// @export
    pub fn d4_track(&self) -> &str {
        self.d4_track.as_deref().unwrap_or("")
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq)]
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
// TODO: returning an R `list!` gives me a mem error and I'm not sure why
//       this is a workaround for now
#[extendr]
impl Query {
    /// Get the chromosome specified by this query
    /// @export
    pub fn chr(&self) -> &str {
        &self.chr
    }

    /// Get the inclusive "left" position specified by this query
    /// @export
    pub fn left(&self) -> u32 {
        self.left
    }

    /// Get the inclusive "right" position specified by this query
    /// @export
    pub fn right(&self) -> u32 {
        self.right
    }
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rd4;
    impl D4File;
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

    use crate::{D4File, Query};

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
        let tempdir = TempDir::new().unwrap();
        let data = create_d4_file(tempdir.path());

        let file = D4File::new(data.to_str().unwrap());
        assert_eq!(file.list_tracks(), vec![String::from("")]);
        assert_eq!(file.chroms(), vec![String::from("chr1")]);
        let result = file.query("chr1", 12, 22);
        assert_eq!(result.results(), &[0, 0, 0, 200, 0, 0, 0, 0, 100, 0]);
        assert_eq!(result.d4_file(), data.to_string_lossy());
        assert_eq!(result.d4_track(), String::from(""));
        let q = result.query();
        assert_eq!(q, &Query::new(String::from("chr1"), 12, 22));
    }
}
