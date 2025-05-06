# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.1] - 2025-05-06

### Fixed

- N2 + O --> NO + N reaction in air_chemistry.json

### Added

- dt=1e-9, fnum=1e5 solution for air chemistry

## [0.9.0] - 2025-04-29

### Changed

- get_nrho to return a dictionary with species names as keys

## [0.8.0] - 2025-03-07

### Changed

- implemented correct relaxation of gas mixtures with temperture dependent relaxation numbers

## [0.7.1] - 2025-03-06

### Fixed

- check for total mole fraction being out of bounds (0, 1]

## [0.7.0] - 2025-02-03

### Changed

- refactored solver to work with number densities instead of mole fractions

## [Unreleased]

## [0.6.1] - 2025-01-31

### Fixed

- minor bugs in output

## [0.6.0] - 2025-01-30

### Added

- reader for SPARTA log files

## [0.5.0] - 2025-01-29

### Added

- calculation of temperature dependant collision numbers for diatomic species

## [0.4.0] - 2025-01-23

### Added

- csv writer
- netCDF4 writer

## [0.3.1] - 2024-12-21

### Fixed

- bug in interface when adding reactions

## [0.3.0] - 2024-12-21

### Added

- Arrhenius type chemical reactions

## [0.2.0] - 2024-12-17

### Added

- Landau-Teller typ relaxation of mixtures with common rotationa and vibrational temperature

## [0.1.0]

### Added

- base project structure