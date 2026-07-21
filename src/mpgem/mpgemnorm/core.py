import os
import logging
import numpy as np
import pandas as pd
from scipy import stats
from importlib import resources


def _data(name):
    """Absolute path to a file bundled in qnorm/data/."""
    return str(resources.files("mpgem.mpgemnorm").joinpath("data", name))


def _fold(s):
    return str(s).upper()


def _get_logger(log_path):
    lg = logging.getLogger("qnorm"); lg.setLevel(logging.INFO); lg.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S")
    for h in (logging.FileHandler(log_path), logging.StreamHandler()):
        h.setFormatter(fmt); lg.addHandler(h)
    return lg


def _save_list(names, path, header):
    pd.Series(list(names), name=header).to_csv(path, sep="\t", index=False)


# ============================================================ step 1: read input
def load_matrix(matrix_path, genes_path=None, samples_path=None, sep="\t"):
    """TSV -> labelled DataFrame; .npy -> array + gene names (+ optional sample ids)."""
    if str(matrix_path).lower().endswith(".npy"):
        arr = np.load(matrix_path)
        if genes_path is None:
            raise ValueError(".npy input requires genes_path (gene names in column order)")
        genes = pd.read_csv(genes_path, sep=sep, header=None)[0].tolist()
        if arr.shape[1] != len(genes):
            raise ValueError(f"matrix has {arr.shape[1]} columns but {len(genes)} gene names")
        if samples_path is not None:
            samples = pd.read_csv(samples_path, sep=sep, header=None)[0].tolist()
            if arr.shape[0] != len(samples):
                raise ValueError(f"matrix has {arr.shape[0]} rows but {len(samples)} sample ids")
        else:
            samples = [f"S{i}" for i in range(arr.shape[0])]
        return pd.DataFrame(arr, index=samples, columns=genes).astype(np.float32)
    return pd.read_csv(matrix_path, sep=sep, index_col=0).astype(np.float32)


# ============================================================ step 2: reference
def load_map(path, sep="\t"):
    """probe -> gene Series from the first two columns (map TSV has a header row)."""
    m = pd.read_csv(path, sep=sep, dtype=str)
    probe_col, gene_col = m.columns[:2]
    m = m[[probe_col, gene_col]].dropna().drop_duplicates(subset=probe_col)
    return m.set_index(probe_col)[gene_col]


def load_gene_averages(path, sep="\t"):
    """Reference table [gene, average]: the reference gene set (order preserved) and the
    per-gene averages of the reference matrix, used to build a subset RSQD."""
    t = pd.read_csv(path, sep=sep)
    genes = t.iloc[:, 0].tolist()
    avg = t.iloc[:, 1].to_numpy(dtype=np.float32)
    return genes, dict(zip(genes, avg))


def load_rqd(path):
    """Precomputed RQD stored as .npy."""
    return np.load(path).astype(np.float32)


# ============================================================ step 3: detect + collapse
def is_probe_level(columns, folded_probes):
    return any(_fold(c) in folded_probes for c in columns)


def collapse_probes_to_genes(expr, folded_probe_to_gene, out_dir, log):
    """Drop probes not in the table (list saved), then max per gene."""
    kept = [c for c in expr.columns if _fold(c) in folded_probe_to_gene]
    foreign = [c for c in expr.columns if _fold(c) not in folded_probe_to_gene]
    if foreign:
        p = os.path.join(out_dir, "dropped_probes.tsv"); _save_list(foreign, p, "dropped_probe")
        log.warning("%d probe(s) not in the database were dropped; list -> %s", len(foreign), p)
    if not kept:
        log.error("no input probes match the database"); raise ValueError("no probes match the table")
    genes = [folded_probe_to_gene[_fold(c)] for c in kept]
    g = expr[kept].T.groupby(genes, sort=True).max().T.astype(np.float32)
    log.info("collapsed %d probes -> %d genes (max per gene)", len(kept), g.shape[1])
    return g


# ============================================================ step 5: NaN
def handle_nan(user, nan_action, out_dir, log):
    """Handle NaN in the gene-level matrix: drop NaN gene columns (default, list saved),
    or fill NaN with zero (nan_action='zero')."""
    n_nan = int(np.isnan(user.to_numpy()).sum())
    if n_nan == 0:
        log.info("NaN check: none found"); return user
    log.warning("NaN check: %d NaN value(s) in %d gene(s)", n_nan, int(user.isna().any().sum()))
    if nan_action == "zero":
        user = user.fillna(0.0)
        log.info("NaN handling: filled NaN with zero")
    else:  # default: drop the NaN-containing genes (list saved)
        bad = user.columns[user.isna().any()].tolist()
        p = os.path.join(out_dir, "nan_dropped_genes.tsv"); _save_list(bad, p, "nan_dropped_gene")
        user = user.drop(columns=bad)
        log.info("NaN handling: dropped %d gene(s) with NaN; list -> %s", len(bad), p)
    return user


# ============================================================ step 7: apply
def assign_quantile(data, quantile, requantile):
    for i in range(0, data.shape[0]):
        requantile[i] = quantile[(stats.rankdata(data[i]) - 1).astype('int')]


# ============================================================ pipeline
def run(matrix_path="user_matrix.tsv",
        rqd_path=None,          # full-overlap runs      (default: bundled)
        ref_avg_path=None,      # gene set + subset RSQD  (default: bundled)
        map_path=None,          # probe -> gene           (default: bundled)
        out_path="normalized_matrix.tsv",
        genes_path=None, samples_path=None,          # for .npy input
        nan_action="drop",
        log_path=None, sep="\t", chunk_rows=200):
    if rqd_path is None:      rqd_path     = _data("reference_rqd.npy")
    if ref_avg_path is None:  ref_avg_path = _data("reference_gene_averages.tsv")
    if map_path is None:      map_path     = _data("probe_map.tsv")

    out_dir = os.path.dirname(out_path) or "."
    os.makedirs(out_dir, exist_ok=True)
    if log_path is None:
        log_path = os.path.join(out_dir, "qnorm.log")
    log = _get_logger(log_path); log.info("=== qnorm start ===")

    # 1. read input
    user = load_matrix(matrix_path, genes_path, samples_path, sep)
    log.info("loaded input: %d samples x %d columns from %s", user.shape[0], user.shape[1], matrix_path)

    # 2. reference vocabulary
    pgmap = load_map(map_path, sep)
    folded_probe_to_gene = {_fold(p): g for p, g in pgmap.items()}
    folded_probes = set(folded_probe_to_gene)
    ref_genes, gene_avg = load_gene_averages(ref_avg_path, sep)
    folded_gene_to_canon = {_fold(g): g for g in ref_genes}
    folded_genes = set(folded_gene_to_canon)
    log.info("reference: %d probes, %d genes (gene averages loaded)", len(folded_probes), len(ref_genes))

    # 3. detect + collapse if probe-level
    if is_probe_level(user.columns, folded_probes):
        log.info("detected: PROBE-level input")
        user = collapse_probes_to_genes(user, folded_probe_to_gene, out_dir, log)
    else:
        log.info("detected: GENE-level input")

    # 4. drop genes not in the reference
    keep_cols = [c for c in user.columns if _fold(c) in folded_genes]
    foreign = [c for c in user.columns if _fold(c) not in folded_genes]
    if foreign:
        p = os.path.join(out_dir, "dropped_genes.tsv"); _save_list(foreign, p, "dropped_gene")
        log.warning("%d gene(s) not in the database were dropped; list -> %s", len(foreign), p)
    if not keep_cols:
        log.error("no genes shared with the database"); raise ValueError("no genes shared with the table")
    user = user[keep_cols]
    log.info("kept %d gene(s) present in the database", len(keep_cols))

    # 5. NaN
    user = handle_nan(user, nan_action, out_dir, log)
    keep_cols = list(user.columns)
    canon = [folded_gene_to_canon[_fold(c)] for c in keep_cols]

    # 6. choose the distribution
    n_common = len(set(canon))
    if n_common == len(ref_genes):
        dist = load_rqd(rqd_path)
        if len(dist) != len(keep_cols):
            log.error("RQD length %d != gene count %d", len(dist), len(keep_cols))
            raise ValueError("RQD length != reference gene count")
        log.info("coverage: FULL -> precomputed RQD (%d genes)", len(keep_cols))
    else:
        dist = np.sort(np.array([gene_avg[g] for g in canon], dtype=np.float32))  # subset RSQD
        log.info("coverage: SUBSET (%d of %d) -> RSQD from per-gene averages", n_common, len(ref_genes))

    # 7. apply + stream output
    first = True
    for start in range(0, user.shape[0], chunk_rows):
        block = user.iloc[start:start + chunk_rows].to_numpy(dtype=np.float32)
        requant = np.empty_like(block)
        assign_quantile(block, dist, requant)
        pd.DataFrame(requant, index=user.index[start:start + chunk_rows], columns=keep_cols).to_csv(
            out_path, sep=sep, mode="w" if first else "a", header=first)
        first = False

    log.info("wrote normalized matrix: %d samples x %d genes -> %s", user.shape[0], len(keep_cols), out_path)
    log.info("=== qnorm done ==="); return out_path


if __name__ == "__main__":
    run()
