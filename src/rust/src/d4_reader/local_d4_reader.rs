//! Implementation of [`D4Reader`] for files on a filesystem.
use std::any::Any;

use d4::stab::SecondaryTablePartReader;
use d4::task::{Histogram, Task, TaskContext};
use d4::Chrom;
use d4::{ptab::DecodeResult, task::Mean};
use ordered_float::OrderedFloat;

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
                                st_value as f64
                            } else {
                                value as f64
                            }
                        }
                    })
            })
            .collect();

        let result = if denominator == 1.0 {
            result.into_iter().map(|v| v as f64).collect()
        } else {
            result.into_iter().map(|v| v as f64 / denominator).collect()
        };

        QueryResult::new(
            query.clone(),
            self.path.clone(),
            track.map(|x| x.to_owned()),
            result,
            None,
        )
    }

    /// Compute the mean depth for each region in `regions` for the specified `track`.
    fn mean(&self, regions: &[Query], track: Option<&str>) -> Vec<f64> {
        let mut reader = self.open(track);
        let result = Mean::create_task(
            &mut reader,
            &regions.iter().map(|r| r.into()).collect::<Vec<(&str, u32, u32)>>(),
        )
        .expect("Failed to run mean tasks")
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
    fn adjust_bin_size(&self, bin_size: u32, allow_bin_size_adjustment: bool) -> u32 {
        bin_size
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

    /// Compute the percentile value for each input region
    pub fn percentile(&self, regions: &[Query], track: Option<&str>, percentile: f64) -> Vec<f64> {
        let results = self.histogram(regions, 0, 1000, track);
        let mut output = vec![];
        for (result, query) in results.iter().zip(regions.iter()) {
            // Inner loop is _percentile_impl in the python code
            // TODO, what is this doing?
            // TODO extract out consts
            // TODO: take ref to query? cloneing the string is enefficient
            let (hist, below, above) = &self.histogram(&[query.clone()], 0, 65536, track)[0];
            let above = *above as f64;
            let below = *below as f64;
            let (chr, begin, end) = query.into();
            let total = (end - begin) as f64;
            if (percentile < below * 100.0 / total) || (100.0 - above * 100.0 / total < percentile)
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
                let mut acc = below;
                for (value, count) in hist {
                    if (acc + *count as f64) * 100.0 / total > percentile {
                        output.push(*value as f64);
                        break;
                    }
                    acc += *count as f64;
                }
            }
        }
        output
    }

    /// histogram(regions, min, max)
    /// --
    ///
    /// Returns the histogram of values in the given regions
    ///
    /// regions: The list of regions we are asking
    /// min: The first bucket of the histogram
    /// max: The exclusive last bucket of this histogram (i.e. a max of 100 will mean the histogram goes up to 99)
    ///
    /// The return value is a list of histograms (including the count of below min and above max
    /// items)
    pub fn histogram(
        &self,
        regions: &[Query],
        min: i32,
        max: i32,
        track: Option<&str>,
    ) -> Vec<(Vec<(i32, u32)>, u32, u32)> {
        // TODO: have this return a Histogram object
        let mut reader = self.open(track);
        let spec = regions
            .iter()
            .map(|r| r.into())
            .map(|(chr, start, stop)| Histogram::with_bin_range(chr, start, stop, min..max))
            .collect();
        let result =
            TaskContext::new(&mut reader, spec).expect("Failed to compute histogram").run();
        let mut output = vec![];
        for item in &result {
            eprintln!("Result Item: {:?}, {:?}, {:?}", item.output.0, item.output.1, item.output.2);
            // The result item is the (number of values less than min, an array with count for each value from min..max, then the number of values >= max)
            let (below, hist, above) = item.output;
            let hist: Vec<_> = hist.iter().enumerate().map(|(a, &b)| (a as i32 + min, b)).collect(); // TODO: I _think_ that min needs to be added to `a` here (this might be a bug in pyd4)
            output.push((hist, *below, *above));
        }
        // Vec<(Vec<(value, count)>, below, above)>
        output
    }
}
