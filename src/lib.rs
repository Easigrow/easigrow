//! fatigue
//! A library of useful fatigue functions
//! Paul White (2015)

pub static COMMENT: &str = "#  ";

extern crate svg;
extern crate log;

#[macro_use]
extern crate lazy_static;

pub mod beta;
pub mod cycle;
pub mod dadn;
pub mod grow;
pub mod tag;
pub mod io;
pub mod plastic;
pub mod material;
pub mod fracto;
pub mod numbers;
pub mod table;
