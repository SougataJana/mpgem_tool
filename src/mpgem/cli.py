"""`mpgem <tool> ...` — dispatches to each tool's subcommand."""
import argparse
from . import __version__
from .mpgemnorm import cli as norm_cli


def main():
    p = argparse.ArgumentParser(prog="mpgem",
                                description="Tools for gene expression matrix processing.")
    p.add_argument("--version", action="version", version=f"mpgem {__version__}")
    sub = p.add_subparsers(dest="tool", required=True, metavar="<tool>")

    qn = sub.add_parser("qnorm", help="reference quantile normalization (GPL570)")
    norm_cli.add_arguments(qn)
    qn.set_defaults(func=norm_cli.main)
    # future tools register here the same way

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
