"""`mpgem qnorm ...` subcommand."""
from .core import run


def add_arguments(p):
    p.add_argument("--matrix", required=True, help="input matrix (TSV, or .npy)")
    p.add_argument("--out", required=True, help="output TSV path")
    p.add_argument("--map", default=None, help="probe->gene TSV (default: bundled)")
    p.add_argument("--rqd", default=None, help="full RQD .npy (default: bundled)")
    p.add_argument("--averages", default=None, help="gene-averages TSV (default: bundled)")
    p.add_argument("--gene-names", default=None, help="gene names file, for .npy input")
    p.add_argument("--sample-ids", default=None, help="sample ids file, for .npy input")
    p.add_argument("--nan-action", default="drop", choices=["drop", "zero"],
                   help="drop NaN genes (default) or fill with zero")
    p.add_argument("--chunk-rows", type=int, default=200)
    return p


def main(a):
    run(matrix_path=a.matrix, out_path=a.out,
        map_path=a.map, rqd_path=a.rqd, ref_avg_path=a.averages,
        genes_path=a.gene_names, samples_path=a.sample_ids,
        nan_action=a.nan_action, chunk_rows=a.chunk_rows)
