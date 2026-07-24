# mpgem

_This repository is a work in progress_ 

Tools for microarray / gene expression matrix processing.

The package currently provides one tool, **mpgemnorm**: reference quantile
normalization of an expression matrix onto a fixed GPL570 reference distribution.
Normalizing to a shared external reference (rather than within each dataset) makes
values comparable across studies. All reference data is bundled with the package —
nothing to download.

On the command line the tool is invoked as `mpgem qnorm`; in Python it is imported as
`mpgem.mpgemnorm`.

## Install

From GitHub:

```
pip install "git+https://github.com/shandarcodes/MPGEM_NORM.git"
```

numpy, pandas and scipy are installed automatically. Requires Python >= 3.10.
This creates the `mpgem` command and bundles the reference data inside the package.

Optionally, in an isolated environment:

```
conda create -n mpgem python=3.12 -y
conda activate mpgem
pip install "git+https://github.com/shandarcodes/MPGEM_NORM.git"
```

Check the install:

```
mpgem --version
mpgem qnorm --help
```

## Usage

### Command line

```
mpgem qnorm --matrix mydata.tsv --out result.tsv
```

Only `--matrix` and `--out` are required; the probe map, RQD, and gene-averages table
are read from the bundled package data. Output and a `qnorm.log` are written to the
folder given in `--out` (the folder is created if it does not exist).

### Python

```python
from mpgem.mpgemnorm import run

run(matrix_path="mydata.tsv", out_path="result.tsv")
```

## Input format

A tab-separated matrix with samples as rows and gene symbols (or GPL570 probe IDs) as
columns; the first column is the sample ID. Columns must be all genes or all probes,
not a mix.

- **Gene-level input** - kept genes are matched against the reference; genes not in the
  reference are dropped (listed in `dropped_genes.tsv`).
- **Probe-level input** - detected automatically, collapsed to gene level (max per gene);
  probes not in the map are dropped (listed in `dropped_probes.tsv`).

NumPy input is also accepted: `--matrix data.npy --gene-names genes.tsv`
(`--sample-ids` optional; defaults to `S0, S1, ...`).

## How it works

For each sample, genes are ranked and each value is replaced by the reference value at
that rank. Two paths are chosen automatically:

| Coverage | Condition | Reference distribution |
|---|---|---|
| FULL | input covers all 19,320 reference genes | precomputed RQD (`reference_rqd.npy`) |
| SUBSET | input covers fewer genes | RSQD built from per-gene averages of those genes |

For the subset case the reference sub-distribution (RSQD) is built by taking the
per-gene averages of the covered genes and sorting them, rather than streaming the full
reference matrix. This holds because gene ranks are stable across the already
quantile-normalized reference.

## Options

```
mpgem qnorm --matrix M.tsv --out DIR/O.tsv                    # gene- or probe-level TSV
mpgem qnorm --matrix M.npy --out DIR/O.tsv --gene-names G.tsv [--sample-ids S.tsv]
mpgem qnorm --matrix M.tsv --out DIR/O.tsv --nan-action zero  # default: drop
mpgem qnorm --matrix M.tsv --out DIR/O.tsv --map ... --rqd ... --averages ...
mpgem qnorm --matrix M.tsv --out DIR/O.tsv --chunk-rows 100   # lower = less memory
```

By default, genes containing NaN are dropped (listed in `nan_dropped_genes.tsv`);
`--nan-action zero` fills them with zero instead. The bundled reference can be
overridden with `--map`, `--rqd`, and `--averages` - supply all three together,
generated from the same reference matrix.

## Reference data

Bundled in `src/mpgem/mpgemnorm/data/`:

- `probe_map.tsv` - GPL570 probe -> gene-symbol map (41,115 probes -> 19,320 genes),
  derived from the GPL570 (Affymetrix Human Genome U133 Plus 2.0) platform annotation.
- `reference_gene_averages.tsv` - per-gene averages of the reference matrix
  (19,320 genes), used to build the subset RSQD.
- `reference_rqd.npy` - the full reference quantile distribution (19,320 values).

The reference statistics were computed from GPL570 expression data obtained from NCBI
GEO. <!-- TODO: add GEO accession(s) and a one-line description of how the reference
matrix was produced. -->

## Citation

<!-- TODO: add once the paper / preprint is available. -->

## License

<!-- TODO: add a LICENSE file (e.g. MIT) and state it here. -->
