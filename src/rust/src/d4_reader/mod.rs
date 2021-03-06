use std::{any::Any, fmt::Debug};

use d4::Chrom;

use crate::{QueryEnv, QueryResultEnv};

pub(crate) mod http_d4_reader;
pub(crate) mod local_d4_reader;

/// 2^16, a constant to use as a marker for high depth.
const HIGH_DEPTH: u32 = 65536;

// `open` can't be defined on this trait since the uderlying D4 reader types
// for HTTP vs Local don't implement a common trait and are two separate types.
// If an associated type were added to this trait to handle this, then the `D4Source`
// would in turn need to be generic over that associated type, even when the `D4Reader`
// is used as a trait object.
//
// Therefore the implementation of open isn't specified here. `open` should only be called
// by these trait functions internally anyways.

/// The trait that defines how a D4 source can be accessed
pub(crate) trait D4Reader: Debug {
    /// Create a [`D4Reader`] that can query the given tracks.
    fn new(source: String) -> Self
    where
        Self: Sized;

    /// Extract the chromosomes from the D4Reader.
    fn get_chroms(&self) -> Vec<Chrom>;

    /// Extract the available tracks in this D4Reader.
    fn get_tracks(&self) -> Vec<String>;

    /// Get the denominator from the header.
    fn get_demoninator(&self) -> f64;

    /// Return the source location for this [`D4Reader`]
    fn source(&self) -> String;

    /// Adjust the bin size selection.
    fn adjust_bin_size(&self, bin_size: u32, allow_bin_size_adjustment: bool) -> u32;

    /// Query the source for the depths over the specified region
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if None, the first track in the source will be queried
    ///
    /// **Note**: This function will change to support queries over multiple tracks at a future time.
    fn query(&self, chr: String, left: u32, right: u32, track: Option<String>) -> QueryResultEnv {
        let query = QueryEnv::new(chr, left, right);
        self.query_track(&query, track.as_deref())
    }

    /// Query a single track with the given [`QueryEnv`].
    ///
    /// This method is intended to be used internally by the `query` method when that method
    /// supports queries over multiple tracks at once.
    fn query_track(&self, query: &QueryEnv, track: Option<&str>) -> QueryResultEnv;

    /// Compute the mean depth for each region in `regions` for the specified `track`.
    fn mean(&self, regions: &[QueryEnv], track: Option<&str>) -> Vec<f64>;

    /// Cast the object implementing this trait to [`Any`] to allow for downcasting.
    fn as_any(&self) -> &dyn Any;
}
