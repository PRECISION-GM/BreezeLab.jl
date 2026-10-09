# Input data

This directory is not versioned (see `.gitignore`). Populate it with

```
julia --project data_wrangling/fetch_covert_inputs.jl
```

which downloads pinned revisions of the *public* SAM configuration files of Covert, Mechem & Zhang (2022,
ACP, doi:10.5194/acp-22-1159-2022) for 18 July 2017 from
<https://github.com/dmechem/ENA_variability_LES_bulk_paper> (bulk paper) and
<https://github.com/dmechem/ENA_variability_LES_bin_paper> (bin paper, with `prm`), and
records their commit hashes and SHA-256 checksums in `data/MANIFEST.txt`.

## Official LASSO-ENA `samin` bundle (ARM account required)

The official experiment definition is the `samin` tar archive of the selected ensemble
member (candidate: `20170718era5s1n0d25x100_sbmwrm-aer2-flxsst`), obtained with a free ARM
account from the LASSO-ENA Bundle Browser
(<https://lasso-ena.svcs.arm.gov/latest/bundle_browser.html>). **Do not commit it.**
Extract it as `data/lasso/<run-id>/` so that `snd`, `lsf`, `sfc`, `prm` (and `grd`) sit in
that directory, then build the member with `ena_lasso(; member, epoch, bundle_dir)`.
See [the protocol instructions](../docs/cases/ena.md#lasso-ena). Record the run id, DOI
(10.5439/2572661), download date and `sha256sum` of the archive in your own provenance
file. `cases/cli/run_case.jl` records extracted input-file checksums automatically;
add the original archive checksum, DOI, and run ID to your experiment record.

### Direct ARM API download

For archives online in ARM's archive, the [ARM Live Data API](https://armlive.svcs.arm.gov/)
can download an exact filename without a browser order. Load `ARM_USERNAME` and
`ARM_TOKEN` from your existing local secret configuration, then run:

```sh
julia data_wrangling/fetch_arm_inputs.jl \
  enalasso_samin_20170718era5s1n0d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar \
  data/lasso/archives
```

The standalone Julia script requires no model initialization. It checks the tar,
rejects redirects and overwrites, and writes a `.tar.toml` sidecar containing the
SHA-256 checksum, size, retrieval time, DOI, and unauthenticated service address.
It does not extract the archive. Keep credentials outside the repository and avoid
putting token literals in shell history. Existing archives are never overwritten.

Availability is per file. On 6 October 2026, the example above returned HTTP 404
from both ARM download endpoint variants; the bundle browser also rejected its
order. A successful ENA aircraft-data query confirmed that the credentials and
API worked. This is an unresolved input-availability issue, not permission to
substitute Covert inputs or to claim official LASSO reproduction. Files not online
may require ARM's normal ordering/staging process.
