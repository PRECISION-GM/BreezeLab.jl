# Data wrangling

- `fetch_covert_inputs.jl`: pinned public SAM input files for the ENA case.
- `fetch_manifest.jl`: checksum-pinned SEA STARR drivers and TRACER protocol references.
- `fetch_arm_inputs.jl`: one exact ARM LASSO `samin` archive, with credentials kept outside
  Git and error messages. These archives are model inputs, not station observations.

These tools populate ignored `data/` directories. They do not build or run simulations.

## NumericalEarth observations interface

Use NumericalEarth's `Metadata`/`Metadatum` and dataset dispatch for observations
where supported. New dataset support can live **here in BreezeLab**, without a fork:
define a BreezeLab-owned dataset type and extend NumericalEarth's hooks in a Julia
package extension (`ext/BreezeLabNumericalEarthExt.jl`, activated by a weak dependency).
The inspected NumericalEarth revision is
[1a06f3b](https://github.com/NumericalEarth/NumericalEarth.jl/tree/1a06f3bc0606608f4e7a93a8729587204df7e20c).
It has the metadata interface and dataset modules, but no ARM observation adapter.

For the first ARM adapter (MWRRET LWP and VDIS rain), the concrete work is:

1. Own the product/site/facility identity in BreezeLab dataset types. Keep requested
   variable, date window and region in NumericalEarth metadata, not the dataset type.
2. Extend `available_variables`, `dataset_variable_name`, `default_region`,
   `default_download_directory`, `metadata_filename` and unit-conversion hooks.
   Use authenticated ARM discovery for actual filenames: daily files do not always
   start at midnight. Do not invent coverage or file timestamps from a date range.
3. Extend `Downloads.download` for metadata parameterized by these owned dataset
   types, using a shared authenticated ARM transport with no tokens in metadata URLs.
4. Add a station-series reader that retains native UTC samples, retrieval QC,
   uncertainty and DQR flags. Station observations are not automatically regular
   grids suitable for `FieldTimeSeries`; sampling against LES needs an explicit operator.
5. Test discovery, units, QC/missing data and metadata dispatch with offline fixtures;
   keep credential-dependent acquisition out of ordinary CI.

This is the selected extension design, **not an implemented ARM adapter**. No new
one-off observation downloader or NumericalEarth dependency is introduced by this
directory reorganization. The existing LASSO archive downloader remains appropriate
for acquiring model setup bundles; it is not a substitute for the observation adapter.
