//! Implementation of [`D4Reader`] for files on a filesystem.
use d4::ptab::DecodeResult;
use d4::stab::SecondaryTablePartReader;
use d4::Chrom;

use crate::{d4_reader::D4Reader, Query, QueryResult};

#[derive(Clone, Debug)]
pub(crate) struct LocalD4Reader {
    path: String,
}

impl D4Reader for LocalD4Reader {
    fn new(source: String) -> Self
    where
        Self: Sized,
    {
        Self { path: source }
    }

    fn get_chroms(&self) -> Vec<Chrom> {
        self.open(None).header().chrom_list().to_vec()
    }

    fn get_tracks(&self) -> Vec<String> {
        let mut tracks = vec![];
        d4::find_tracks_in_file(&self.path, |_| true, &mut tracks)
            .expect("Failed to extract tracks from file.");
        tracks.into_iter().map(|x| x.to_string_lossy().to_string()).collect()
    }

    fn source(&self) -> String {
        self.path.clone()
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
        let partition =
            reader.split(None).unwrap_or_else(|_| panic!("Failed to partition {:?}", self.path));
        let chr = query.chr.clone();
        let result: Vec<i32> = partition
            .into_iter()
            .flat_map(move |(mut ptab, mut stab)| {
                let (part_chr, begin, end) = ptab.region();
                let part_chr = part_chr.to_string();
                let pd = ptab.to_codec();
                (if part_chr == chr { query.left.max(begin)..query.right.min(end) } else { 0..0 })
                    .map(move |pos| match pd.decode(pos as usize) {
                        DecodeResult::Definitely(value) => value,
                        DecodeResult::Maybe(value) => {
                            if let Some(st_value) = stab.decode(pos) {
                                st_value
                            } else {
                                value
                            }
                        }
                    })
            })
            .collect();
        QueryResult::new(query.clone(), self.path.clone(), track.map(|x| x.to_owned()), result)
    }
}

impl LocalD4Reader {
    /// Open a [`LocalD4Reader`] for the specified track
    fn open(&self, track: Option<&str>) -> d4::D4TrackReader {
        // TODO: D4TrackReader::open splits the path on `:` and takes anything to the right as a track specifier
        // There should be a way to create a reader like remote with the optional track name
        let track_spec =
            if let Some(track) = track { [&self.path, track].join(":") } else { self.path.clone() };
        d4::D4TrackReader::open(&track_spec)
            .unwrap_or_else(|_| panic!("Failed to create local reader for: {:?}", self.path))
    }
}
