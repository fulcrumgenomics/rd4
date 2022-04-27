//! Implementation of [`D4Reader`] for files on a filesystem.
use std::any::Any;

use d4::stab::SecondaryTablePartReader;
use d4::task::{Task, TaskContext};
use d4::Chrom;
use d4::{ptab::DecodeResult, task::Mean};
use ordered_float::OrderedFloat;

use crate::Histogram;
use crate::{d4_reader::D4Reader, Query, QueryResult};

use super::HIGH_DEPTH;

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

    fn query_track(&self, query: &Query, track: Option<&str>) -> QueryResult {
        let mut reader = self.open(track);
        let partition =
            reader.split(None).unwrap_or_else(|_| panic!("Failed to partition {:?}", self.path));
        let chr = query.chr.clone();
        let denominator = self.get_demoninator();
        let result: Vec<f64> = partition
            .into_iter()
            .flat_map(move |(mut ptab, mut stab)| {
                let (part_chr, begin, end) = ptab.region();
                let part_chr = part_chr.to_string();
                let pd = ptab.to_codec();
                (if part_chr == chr { query.left.max(begin)..query.right.min(end) } else { 0..0 })
                    .map(move |pos| match pd.decode(pos as usize) {
                        DecodeResult::Definitely(value) => value as f64,
                        DecodeResult::Maybe(value) => {
                            if let Some(st_value) = stab.decode(pos) {
                                st_value as f64 / denominator
                            } else {
                                value as f64 / denominator
                            }
                        }
                    })
            })
            .collect();

        QueryResult::new(
            query.clone(),
            self.path.clone(),
            track.map(|x| x.to_owned()),
            result,
            None,
        )
    }

    fn mean(&self, regions: &[Query], track: Option<&str>) -> Vec<f64> {
        let mut reader = self.open(track);
        let result = Mean::create_task(
            &mut reader,
            &regions.iter().map(|r| r.into()).collect::<Vec<(&str, u32, u32)>>(),
        )
        .expect("Failed to create mean tasks")
        .run();
        let mut output = vec![];
        for item in &result {
            output.push(*item.output);
        }
        output
    }

    fn percentile(&self, regions: &[Query], track: Option<&str>, percentile: f64) -> Vec<f64> {
        let mut output = vec![];
        // Inner loop is _percentile_impl in the python code
        for query in regions {
            let hist = &self.histogram(&[query.clone()], 0, HIGH_DEPTH as i32, track)[0];
            let above = hist.above as f64;
            let below = hist.below as f64;
            let (_chr, begin, end) = query.into();
            let total = (end - begin) as f64;
            if (percentile < (below / total) * 100.0)
                || (100.0 - ((above / total) * 100.0) < percentile)
            {
                let data = self.query_track(query, track);
                let mut low = data
                    .results()
                    .iter()
                    .copied()
                    .map(OrderedFloat)
                    .min()
                    .unwrap_or(OrderedFloat(0.0))
                    .into_inner();
                let mut high = data
                    .results()
                    .iter()
                    .copied()
                    .map(OrderedFloat)
                    .max()
                    .unwrap_or(OrderedFloat(0.0))
                    .into_inner()
                    + 1.0;
                while high - low > 1.0 {
                    let mid = (high + low) % 2.0;
                    let p = data.results().iter().filter(|depth| **depth < mid).sum::<f64>()
                        * 100.0
                        / total;
                    if p < percentile {
                        low = mid;
                    } else {
                        high = mid;
                    }
                }
                output.push(low);
            } else {
                let mut value = hist.first_value;
                for prefix_sum in hist.prefix_sum.iter().skip(1) {
                    if (*prefix_sum as f64 / total as f64) * 100.0 > percentile {
                        output.push(value as f64);
                        break;
                    }
                    value += 1;
                }
            }
        }
        output
    }

    fn histogram(
        &self,
        regions: &[Query],
        min: i32,
        max: i32,
        track: Option<&str>,
    ) -> Vec<Histogram> {
        let mut reader = self.open(track);
        let spec = regions
            .iter()
            .map(|r| r.into())
            .map(|(chr, start, stop)| {
                d4::task::Histogram::with_bin_range(chr, start, stop, min..max)
            })
            .collect();
        let result =
            TaskContext::new(&mut reader, spec).expect("Failed to compute histogram").run();
        let mut output = vec![];
        for item in &result {
            // The result item is the (number of values less than min, an array with count for each value from min..max, then the number of values >= max)
            let (below, hist, above) = item.output;
            let hist: Vec<_> = hist.iter().enumerate().map(|(a, &b)| (a as i32 + min, b)).collect(); // This differs from the pyd4 api, which does not add `a + min` and just lets the "value" be 0
            let hist = Histogram::new(hist, *below, *above);
            output.push(hist);
        }
        output
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn get_demoninator(&self) -> f64 {
        self.open(None).header().get_denominator()
    }

    /// For local files no adjustment is needed.
    fn adjust_bin_size(&self, bin_size: u32, _allow_bin_size_adjustment: bool) -> u32 {
        bin_size
    }
}

impl LocalD4Reader {
    /// Open a [`LocalD4Reader`] for the specified track
    fn open(&self, track: Option<&str>) -> d4::D4TrackReader {
        let track_spec =
            if let Some(track) = track { [&self.path, track].join(":") } else { self.path.clone() };
        d4::D4TrackReader::open(&track_spec)
            .unwrap_or_else(|_| panic!("Failed to create local reader for: {:?}", self.path))
    }
}
