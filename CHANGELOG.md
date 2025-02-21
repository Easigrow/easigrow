# Changelog - Easigrow
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2025-02-21
### Added
- New material for Linear Lower Threshold 5b dadn.


## [2.0.0] - 2024-11-28
### Added
- Support for different crack growth models.
- Willenborg crack growth model.
- 1D beta table.
### Changed
- `--list` output revised and improved.
- Behaviour of CLI option `--astart`.
- Format for specifying file-based beta definitions. (See beta section of `easigrow --list`)
### Fixed
- All betas have been revised and corrected where needed.


## [1.2.0] - 2022-08-02
### Added
- Log crack depth for `da` in `da/dB` during optimisations.
- CLI option `--opt_linear_crack`.
### Changed
- Crack growth calculation during optimisations defaults to log; was previously linear.
### Removed
- GSL support


## [1.1.0] - 2022-07-08
### Added
- Particle Swarm optimisation.
- CLI options `--bound_min`, `--bound_max`, `--swarm_size` and `--unbounded`.
- DaDn equations: Linear Lower Threshold 5A, 6 and 6A.
- Crate dependencies: rand 0.8.4, cubic-splines 0.2.0.
### Changed
- Optimisation now uses "nearest block" logic.
- `grow::History::grow_crack()` now accepts a reference, resulting in improved performance.
- Refactored `optimisation` into its own directory.
- Refactored `table` into its own directory.
### Fixed
- Improved the optimisation implementation so it no longer crashes.
- Error during optimisation if attempting to load cycles using --cycle_infile.
- Non-GSL table interpolation.
- GSL implementation.
- Various clippy warnings.


## [1.0.0] - 2022-01-13
### Added
- This changelog file.
- DaDn equation: Linear Lower Threshold 5.
- Common input parameters `rmax` and `kneg` (currently only used by the above equation).
- Crate dependency: lazy_static 1.4.0.
- `make_model_from_file()` function in the `dadn` module.
- `Options` struct in the `dadn` module for passing fixed equation parameters.

### Changed
- Command line parameter input format.
- Rust edition to 2018.

### Removed
- `update()` and `variables()` methods from the `DaDn` trait.

### Fixed
- Various clippy warnings about syntax.
- Various code improvements for quality and readability. 
