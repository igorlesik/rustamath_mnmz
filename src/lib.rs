//! Mimimization functions.
//!
//! (c) 2023 Igor Lesik
//! MIT license
//!
//! Task of minimization: for given function _f_ that depends on one or more independent
//! variables, find the value of those variables where _f_ takes on a minimum value.
//!
pub mod bracket;
pub use bracket::{find_bracket, BracketRes};
pub mod golden_section;
pub use golden_section::golden_section_search;
pub mod brents_method;
pub use brents_method::brent_search;
pub mod brents_df_method;
pub use brents_df_method::brent_df_search;
pub mod simplex;
pub use simplex::amoeba;

#[cfg(test)]
#[macro_use]
extern crate assert_float_eq;