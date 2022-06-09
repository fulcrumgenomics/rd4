//! Implementation of [`D4Reader`] for files on a filesystem.
use std::any::Any;
use std::panic;

use d4::stab::SecondaryTablePartReader;
use d4::task::{Task, TaskContext};
use d4::Chrom;
use d4::{ptab::DecodeResult, task::Mean};
use ordered_float::OrderedFloat;

use crate::HistogramEnv;
use crate::{d4_reader::D4Reader, QueryEnv, QueryResultEnv};

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

    fn query_track(&self, query: &QueryEnv, track: Option<&str>) -> QueryResultEnv {
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

        QueryResultEnv::new(
            query.clone(),
            self.path.clone(),
            track.map(|x| x.to_owned()),
            result,
            None,
        )
    }

    fn mean(&self, regions: &[QueryEnv], track: Option<&str>) -> Vec<f64> {
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
        // TODO: D4TrackReader::open splits the path on `:` and takes anything to the right as a track specifier
        // There should be a way to create a reader like remote with the optional track name
        
        let try_path = d4::D4TrackReader::open(&self.path)
            .unwrap_or_else(|_| panic!("Could not create local reader for path: {:?}", self.path));
        
        if let Some(track) = track {
            let track_spec = [&self.path, track].join(":");
            d4::D4TrackReader::open(&track_spec)
                .unwrap_or_else(|_| panic!("Failed to open track {:?} for path: {:?}", track, self.path))
        } else {
            try_path
        }
    }

    /// Compute the percentile value for each input region
    pub fn percentile(&self, regions: &[QueryEnv], track: Option<&str>, percentile: f64) -> Vec<f64> {
        let mut output = vec![];
        // Inner loop is _percentile_impl in the python code
        for query in regions {
            if percentile <= 0.0 {
                output.push(0.0)
            } else {
                let data = self.query_track(query, track);
                let min = data
                    .results()
                    .iter()
                    .copied()
                    .map(OrderedFloat)
                    .min()
                    .unwrap_or(OrderedFloat(0.0))
                    .into_inner();
                let max = data
                    .results()
                    .iter()
                    .copied()
                    .map(OrderedFloat)
                    .max()
                    .unwrap_or(OrderedFloat(0.0))
                    .into_inner();
                if min == max || percentile >= 100.0 {
                    output.push(max);
                } else {
                    let hist = &self.histogram(&[query.clone()], 0, HIGH_DEPTH as i32, track)[0];
                    let total = hist.total_count() as f64;
                    let above = hist.above as f64;
                    let below = hist.below as f64;
                    if (percentile < (below / total) * 100.0)
                        || (100.0 - ((above / total) * 100.0) < percentile) // Percentile lies outside the histogram
                    {
                        let mut low = min.clone();
                        let mut high = max.clone() + 1.0;
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
                    } else { // Percentile lies inside the histogram
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
            }
        }
        output
    }

    /// histogram(regions, min, max)
    ///
    /// Returns the histogram of values in the given regions
    ///
    /// # Arguments
    /// - `regions` - The list of regions we are asking
    /// - `min` - The first bucket of the histogram
    /// - `max` - The exclusive last bucket of this histogram (i.e. a max of 100 will mean the histogram goes up to 99)
    pub fn histogram(
        &self,
        regions: &[QueryEnv],
        min: i32,
        max: i32,
        track: Option<&str>,
    ) -> Vec<HistogramEnv> {
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
            panic!("\nhist first counts: 0:{:?}, 1:{:?}, 2:{:?}, 3:{:?}", hist[0], hist[1], hist[2], hist[3]);
            let hist: Vec<_> = hist.iter().enumerate().map(|(a, &b)| (a as i32 + min, b)).collect(); // TODO - This differs from the pyd4 api, which does not add `a + min` and just lets the "value" be 0
            let hist = HistogramEnv::new(hist, *below, *above);
            output.push(hist);
        }
        output
    }
}
