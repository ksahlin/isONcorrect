//! isONcorrect --- Rust port.
//!
//! The Python implementation under `src/isoncorrect/` is the normative
//! reference: where the two disagree, Python is right until a human decides
//! otherwise. Equivalence is checked, not assumed --- see `bench/` and
//! `PORTING.md`.
//!
//! Current state: **CLI parity only.** Argument parsing, defaults, validation
//! order and exit codes match the reference; the correction algorithm is not
//! ported yet and the binaries exit non-zero when asked to do real work.

pub mod align;
pub mod anchors;
pub mod cli;
pub mod contexts;
pub mod corrections;
pub mod driver;
pub mod editdist;
pub mod fanout;
pub mod fastq;
pub mod guard;
pub mod minimizers;
pub mod msa;
pub mod params;
pub mod parasail;
pub mod poa;
pub mod profile;
pub mod quality;
pub mod regions;
pub mod simd;
pub mod support;
pub mod validate;
pub mod wis;
