//! Implementation of [`D4Reader`] on remote files over HTTP.
use d4::ssio::http::HttpReader;
use d4::Chrom;

use crate::{d4_reader::D4Reader, Query, QueryResult};

/// Access D4 sources over HTTP
#[derive(Clone, Debug)]
pub struct HttpD4Reader {
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

    fn query(&self, chr: String, left: u32, right: u32, track: Option<String>) -> QueryResult {
        let query = Query::new(chr, left, right);

        if let Some(track) = track {
            self.query_track(&query, Some(&track))
        } else {
            self.query_track(&query, None)
        }
    }
    fn query_track(&self, query: &Query, track: Option<&str>) -> QueryResult {
        let mut reader = self.open(track);
        let result = reader
            .get_view(&query.chr, query.left, query.right)
            .unwrap_or_else(|_| {
                panic!(
                    "Failed to get view for {:?} {}:{}-{}",
                    self.source, &query.chr, query.left, query.right
                )
            })
            .map(|res| if let Ok((_, value)) = res { value } else { 0 })
            .collect();
        QueryResult::new(query.clone(), self.source.clone(), track.map(|x| x.to_owned()), result)
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
