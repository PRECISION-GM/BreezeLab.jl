# Input data

This directory is not versioned (see `.gitignore`). Populate it with

```
julia --project scripts/fetch_inputs.jl
```

which downloads pinned revisions of the *public* SAM configuration files of Covert, Mechem & Zhang (2022,
ACP, doi:10.5194/acp-22-1159-2022) for 18 July 2017 from
<https://github.com/dmechem/ENA_variability_LES_bulk_paper> (bulk paper) and
<https://github.com/dmechem/ENA_variability_LES_bin_paper> (bin paper, with `prm`), and
records their commit hashes and SHA-256 checksums in `data/MANIFEST.txt`.

## Official LASSO-ENA `samin` bundle (restricted)

The official experiment definition is the `samin` tar archive of the selected ensemble
member (candidate: `20170718era5s1n0d25x100_sbmwrm-aer2-flxsst`), obtained with a free ARM
account from the LASSO-ENA Bundle Browser
(<https://lasso-ena.svcs.arm.gov/latest/bundle_browser.html>). **Do not commit it.**
Extract it as `data/lasso/<run-id>/` so that `snd`, `lsf`, `sfc`, `prm` (and `grd`) sit in
that directory, then select `ena_simulation(...; protocol=:lasso_ena_official, epoch=...)`.
See [the protocol instructions](../docs/cases/ena.md#lasso-ena). Record the run id, DOI
(10.5439/2572661), download date and `sha256sum` of the archive in your own provenance
file. `scripts/run_case.jl` records extracted input-file checksums automatically;
add the original archive checksum, DOI, and run ID to your experiment record.
