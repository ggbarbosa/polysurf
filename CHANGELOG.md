# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [v0.1.0] - 2025-09-27

### Added
- First functional version of `polysurf`.
- Implemented the `build()` function for automated construction of polymer–surface systems.
- Modular code structure under `src/`, ready for installation via `pip`.
- Support for:
  - Automatic polymer stretching to match surface lattice dimensions.
  - Vertical alignment with customizable spacing (`v_step`, `z_void`).
  - Automatic orthogonal replication to ensure minimum lateral distance (`lateral_min`).
- Support for `.xyz` file formats via ASE.
- Output files:
  - `slab_super.xyz`
  - `polymer_stretched.xyz`
  - `combined.xyz`
- Metadata export to `metadata.json`.

---
