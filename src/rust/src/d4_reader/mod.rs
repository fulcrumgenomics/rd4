use std::fmt::Debug;

use d4::Chrom;

use crate::{Query, QueryResult};

pub(crate) mod http_d4_reader;
pub(crate) mod local_d4_reader;

// `open` can't be defined on this trait since the uderlying D4 reader types
// for HTTP vs Local don't implement a common trait and are two separate types.
// If an associated type were added to this trait to handle this, then the `D4Source`
// would in turn need to be generic over that associated type, even when the `D4Reader`
// is used as a trait object.
//
// Therefore the implementation of open isn't specified here. `open` should only be called
// by these trait functions internally anyways.

/// The trait that defines how a D4 source can be accessed
pub trait D4Reader: Debug {
    /// Create a [`D4Reader`] that can query the given tracks.
    fn new(source: String) -> Self
    where
        Self: Sized;

    /// Extract the chromosomes from the D4Reader.
    fn get_chroms(&self) -> Vec<Chrom>;

    /// Extract the available tracks in this D4Reader.
    fn get_tracks(&self) -> Vec<String>;

    /// Return the source location for this [`D4Reader`]
    fn source(&self) -> String;

    /// Query the source for the depths over the specified region
    ///
    /// # Arguments
    /// - `chr` - the chromosome to query
    /// - `left` - the inclusive start position
    /// - `right` - the exclusive end position
    /// - `track` - the optional track to query, if None, the first track in the source will be queried
    ///
    /// **Note**: This function will change to support queries over multiple tracks at a future time.
    fn query(&self, chr: String, left: u32, right: u32, track: Option<String>) -> QueryResult;

    /// Query a single track with the given [`Query`].
    ///
    /// This method is intended to be used internally by the `query` method when that method
    /// supports queries over multiple tracks at once.
    fn query_track(&self, query: &Query, track: Option<&str>) -> QueryResult;
}
