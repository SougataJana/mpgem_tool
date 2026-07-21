# qnorm — commands

## 1. Install

```
conda create -n qnorm-test python=3.12 -y
conda activate qnorm-test
cd <folder containing pyproject.toml>
pip install -e .
```

Installs numpy, pandas and scipy automatically and creates the `qnorm` command.

## 2. Check the install

```
which qnorm
qnorm --help
python -c "from qnorm.core import _data; print(_data('probe_map.tsv'))"
```

The last command prints the path of the bundled reference data inside the package.

## 3. Run — gene-level test matrix

`demo_matrix.tsv`: 5 samples x 2,000 genes, gene symbols as columns.

```
mkdir -p ~/qnorm_test
cd <folder containing demo_matrix.tsv>
qnorm --matrix demo_matrix.tsv --out ~/qnorm_test/out.tsv
ls -l ~/qnorm_test/
```

Only `--matrix` and `--out` are needed; the probe map, RQD and gene-averages table come
from the bundled package data.

Expected log:

```
loaded input: 5 samples x 2000 columns from demo_matrix.tsv
reference: 41115 probes, 19320 genes (gene averages loaded)
detected: GENE-level input
kept 2000 gene(s) present in the database
NaN check: none found
coverage: SUBSET (2000 of 19320) -> RSQD from per-gene averages
wrote normalized matrix: 5 samples x 2000 genes -> ~/qnorm_test/out.tsv
```

Produces `out.tsv` and `qnorm.log` in `~/qnorm_test/`.

## 4. Run — probe-level test matrix (with unknown probes)

`demo_probe_matrix.tsv`: 5 samples x 1,020 columns = 1,000 real GPL570 probe IDs plus 20
probe IDs that are not in the database. Exercises probe detection, dropping of unknown
probes, and collapsing probes to genes.

```
mkdir -p ~/qnorm_test_probe
qnorm --matrix demo_probe_matrix.tsv --out ~/qnorm_test_probe/out.tsv
ls -l ~/qnorm_test_probe/
head -3 ~/qnorm_test_probe/dropped_probes.tsv
```

Expected log:

```
loaded input: 5 samples x 1020 columns from demo_probe_matrix.tsv
detected: PROBE-level input
WARNING | 20 probe(s) not in the database were dropped; list -> ~/qnorm_test_probe/dropped_probes.tsv
collapsed 1000 probes -> 973 genes (max per gene)
kept 973 gene(s) present in the database
coverage: SUBSET (973 of 19320) -> RSQD from per-gene averages
wrote normalized matrix: 5 samples x 973 genes -> ~/qnorm_test_probe/out.tsv
```

Produces `out.tsv`, `qnorm.log` and `dropped_probes.tsv` in `~/qnorm_test_probe/`.
The 20 dropped entries are the `FAKE_PROBE_*_at` columns. Output columns are gene
symbols, not probe IDs.

## 5. Where the files are saved

Everything is written to **the folder given in `--out`** — never inside the package.

| File | Example path | When |
|---|---|---|
| normalized matrix | `~/qnorm_test/out.tsv` | always |
| log | `~/qnorm_test/qnorm.log` | always |
| dropped probes | `~/qnorm_test_probe/dropped_probes.tsv` | probe-level input with unknown probes (section 4) |
| dropped genes | `<out folder>/dropped_genes.tsv` | gene-level input with genes not in the reference |
| NaN-dropped genes | `<out folder>/nan_dropped_genes.tsv` | input has NaN and `--nan-action drop` |

Notes:

- The output folder is created automatically if it does not exist.
- The log is written to the file **and** printed to the terminal; both contain the same lines.
- The log is always named `qnorm.log` and is overwritten by the next run that writes to the
  same folder — use a different `--out` folder per run to keep logs (as in sections 3 and 4).
- If `--out out.tsv` is given with no folder, everything lands in the current directory.

Inspect what was produced:

```
ls -l ~/qnorm_test/ ~/qnorm_test_probe/
cat ~/qnorm_test/qnorm.log
```

## 6. Verify the results

```
cat > ~/check_qnorm.py << 'PYEOF'
import os, pandas as pd, numpy as np
from qnorm.core import _data

a = pd.read_csv(_data("reference_gene_averages.tsv"), sep="\t")
ga = dict(zip(a.gene, a.average))

def check(label, path, n_expected):
    o = pd.read_csv(os.path.expanduser(path), sep="\t", index_col=0)
    expected = np.sort(np.array([ga[g] for g in o.columns], dtype=np.float32))
    sorted_rows = np.sort(o.to_numpy(), axis=1)
    print(f"{label}")
    print("  shape:", o.shape, "| genes ==", n_expected, ":", o.shape[1] == n_expected)
    print("  sorts to subset RSQD:", bool(np.allclose(sorted_rows, expected, atol=1e-4)))
    print("  all samples identical when sorted:", bool(np.allclose(sorted_rows.std(axis=0), 0)))

check("gene-level  (section 3)", "~/qnorm_test/out.tsv", 2000)
check("probe-level (section 4)", "~/qnorm_test_probe/out.tsv", 973)

d = pd.read_csv(os.path.expanduser("~/qnorm_test_probe/dropped_probes.tsv"), sep="\t")
print("dropped probes:", len(d), "| all FAKE_:", d.iloc[:, 0].str.startswith("FAKE_PROBE_").all())
PYEOF
python ~/check_qnorm.py
```

Expected:

```
gene-level  (section 3)
  shape: (5, 2000) | genes == 2000 : True
  sorts to subset RSQD: True
  all samples identical when sorted: True
probe-level (section 4)
  shape: (5, 973) | genes == 973 : True
  sorts to subset RSQD: True
  all samples identical when sorted: True
dropped probes: 20 | all FAKE_: True
```

`sorts to subset RSQD: True` means every sample was mapped onto the reference
distribution restricted to its genes. `all samples identical when sorted: True` means the
distributions were equalized across samples.

## 7. Options

```
qnorm --matrix M.tsv --out DIR/O.tsv                    # gene- or probe-level TSV
qnorm --matrix M.npy --out DIR/O.tsv --gene-names G.tsv [--sample-ids S.tsv]
qnorm --matrix M.tsv --out DIR/O.tsv --nan-action zero  # default: drop
qnorm --matrix M.tsv --out DIR/O.tsv --map ... --rqd ... --averages ...
qnorm --matrix M.tsv --out DIR/O.tsv --chunk-rows 100   # lower = less RAM
```

Input format: samples as rows, gene symbols or GPL570 probe IDs as columns, first column
= sample ID. Columns must be all genes or all probes, not a mix. For `.npy` input a
gene-names file is required (one name per line, matching column order); sample IDs are
optional and default to `S0, S1, ...`.

## 8. Uninstall

```
pip uninstall qnorm -y
conda deactivate
conda env remove -n qnorm-test -y
```
