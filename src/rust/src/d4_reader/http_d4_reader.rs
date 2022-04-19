//! Implementation of [`D4Reader`] on remote files over HTTP.
use std::any::Any;

use d4::ssio::http::HttpReader;
use d4::Chrom;

use crate::{d4_reader::D4Reader, Query, QueryResult};

use super::HIGH_DEPTH;

/// Access D4 sources over HTTP
#[derive(Clone, Debug)]
pub(crate) struct HttpD4Reader {
    /// The URL of the D4 source
    source: String,
}

impl D4Reader for HttpD4Reader {
    fn new(source: String) -> Self {
        Self { source }
    }

    fn get_chroms(&self) -> Vec<Chrom> {
        self.open(None).chrom_list().to_vec()
    }

    fn get_tracks(&self) -> Vec<String> {
        let mut tracks = vec![];
        d4::find_tracks_in_file(&self.source, |_| true, &mut tracks)
            .expect("Failed to extract tracks from file.");
        tracks.into_iter().map(|x| x.to_string_lossy().to_string()).collect()
    }

    fn source(&self) -> String {
        self.source.clone()
    }

    fn query_track(&self, query: &Query, track: Option<&str>) -> QueryResult {
        let mut reader = self.open(track);
        let denominator = self.get_demoninator();
        let result: Vec<f64> = reader
            .get_view(&query.chr, query.left, query.right)
            .unwrap_or_else(|_| {
                panic!(
                    "Failed to get view for {:?} {}:{}-{}",
                    self.source, &query.chr, query.left, query.right
                )
            })
            .map(|res| if let Ok((_, value)) = res { value as f64 / denominator } else { 0.0 })
            .collect();

        QueryResult::new(
            query.clone(),
            self.source.clone(),
            track.map(|x| x.to_owned()),
            result,
            None,
        )
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn mean(&self, regions: &[Query], track: Option<&str>) -> Vec<f64> {
        let mut reader = self.open(track);
        let mut output = Vec::with_capacity(regions.len());
        let index =
            reader.load_data_index::<d4::index::Sum>().expect("Failed to load remote index");
        for (chr, begin, end) in regions.iter().map(|r| r.into()) {
            let index_res = index.query(chr, begin, end).expect("Failed to run remote query");
            let sum_res = index_res
                .get_result(&mut reader)
                .expect("Failed to extract query result from remote query.");
            let mean = sum_res.mean(index_res.query_size());
            output.push(mean);
        }
        output
    }

    fn get_demoninator(&self) -> f64 {
        self.open(None).get_header().get_denominator()
    }

    /// Try to adjust the bin size if allowed.
    fn adjust_bin_size(&self, bin_size: u32, allow_bin_size_adjustment: bool) -> u32 {
        if allow_bin_size_adjustment {
            if bin_size < HIGH_DEPTH {
                HIGH_DEPTH
            } else {
                bin_size - bin_size % (HIGH_DEPTH)
            }
        } else {
            bin_size
        }
    }
}

impl HttpD4Reader {
    /// Open an [`HttpD4Reader`] on the specified track.
    fn open(&self, track: Option<&str>) -> d4::ssio::D4TrackReader<HttpReader> {
        // TODO: D4TrackReader::open splits the path on `:` and takes anything to the right as a track specifier
        // There should be a way to create a reader like remote with the optional track name
        let conn = d4::ssio::http::HttpReader::new(&self.source)
            .unwrap_or_else(|_| panic!("Failed to open remote reader for: {:?} ", self.source));
        d4::ssio::D4TrackReader::from_reader(conn, track).unwrap_or_else(|_| {
            panic!("Failed create remote reader for: {:?}:{:?} ", self.source, track)
        })
    }
}
