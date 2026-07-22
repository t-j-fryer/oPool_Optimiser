#!/usr/bin/env python3
"""Friendly command-line interface for the oPool cloning workflow."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from opool_workflow import BUILTIN_CODON_SPECIES, WorkflowConfig, run_workflow


def positive_integer(value: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error
    if result < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def optional_positive_integer(value: str) -> int | None:
    if value.strip().lower() in {"auto", "none", "unlimited"}:
        return None
    return positive_integer(value)


def optional_seed(value: str) -> int | None:
    if value.strip().lower() in {"none", "random"}:
        return None
    try:
        return int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer or 'random'") from error


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="opool",
        description=(
            "Optimize genes, assign fast Golden Gate fragment pools, and add "
            "orthogonal primers/Type IIS elements. Input may be a two-column "
            "amino-acid CSV or an existing *_Optimised.csv file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Example:\n"
            "  python scripts/opool_cli.py --input /path/opTF010_Optimised.csv "
            "--overhangs /path/overhangs_bgal.csv --opool-length 350 "
            "--vector-oh1 TATG --vector-oh2 GGAT --genes-per-subpool 1 --force"
        ),
    )

    essentials = parser.add_argument_group("essential inputs")
    essentials.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Amino-acid CSV or optimized CSV; type is detected automatically.",
    )
    essentials.add_argument(
        "--opool-length",
        type=positive_integer,
        default=250,
        help="Final oligo length in nucleotides.",
    )
    essentials.add_argument(
        "--overhangs",
        type=Path,
        default=None,
        help="Overhang inventory CSV; defaults to data/overhangs.csv.",
    )
    essentials.add_argument(
        "--vector-oh1",
        default="GCTT",
        help="Left vector overhang (4 nt).",
    )
    essentials.add_argument(
        "--vector-oh2",
        default="AGTG",
        help="Right vector overhang (4 nt).",
    )
    essentials.add_argument(
        "--genes-per-subpool",
        type=optional_positive_integer,
        default=None,
        metavar="N|auto",
        help="Maximum genes per sub-pool; use 1 for one gene per block.",
    )

    outputs = parser.add_argument_group("output and input handling")
    outputs.add_argument(
        "--input-kind",
        choices=["auto", "aa", "optimized"],
        default="auto",
        help="Override automatic input detection.",
    )
    outputs.add_argument("--output-dir", type=Path, default=None, help="Output directory; defaults beside input.")
    outputs.add_argument("--run-name", default=None, help="Output filename prefix; defaults to input folder name.")
    outputs.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    outputs.add_argument("--quiet", action="store_true", help="Hide progress messages.")

    assembly = parser.add_argument_group("assembly options")
    assembly.add_argument(
        "--short-pool-max-size",
        type=optional_positive_integer,
        default=1000,
        metavar="N|unlimited",
        help="Separate cap for single-fragment genes.",
    )
    assembly.add_argument("--max-fragments", type=positive_integer, default=32)
    assembly.add_argument(
        "--keep-nterm-met",
        action="store_true",
        help="Keep a leading ATG; by default it is stripped before fragmentation.",
    )
    assembly.add_argument("--typeiis-site", default="GGTCTC", help="Recognition sequence without the N base.")
    assembly.add_argument("--typeiis-n", default="A", help="Single N-position base appended to the Type IIS site.")

    primers = parser.add_argument_group("primer and stuffer options")
    primers.add_argument(
        "--primers",
        type=Path,
        default=None,
        help="Orthogonal-primer CSV; defaults to data/orthogonal_oligos.csv.",
    )
    primers.add_argument(
        "--primer-mode",
        choices=["combinatorial", "unique_pairs"],
        default="combinatorial",
    )
    primers.add_argument("--primer-start-at", default=None, help="First primer ID for unique_pairs mode.")
    primers.add_argument("--skip-primers", action="store_true", help="Stop after pool assignment.")
    primers.add_argument("--no-stuffer", action="store_true", help="Do not pad oligos to opool-length.")
    primers.add_argument(
        "--allow-unsafe-primers",
        action="store_true",
        help="Allow primers containing the Type IIS site or reverse complement.",
    )
    primers.add_argument("--seed", type=optional_seed, default=123, metavar="INT|random")

    optimization = parser.add_argument_group("gene-optimization options (amino-acid inputs only)")
    optimization.add_argument("--gc-min", type=float, default=0.30)
    optimization.add_argument("--gc-max", type=float, default=0.70)
    optimization.add_argument("--gc-window", type=positive_integer, default=50)
    optimization.add_argument(
        "--codon-species",
        default="e_coli",
        metavar="SPECIES|TAXID",
        help=(
            "Codon-optimization host. Built-in keywords: "
            + ", ".join(BUILTIN_CODON_SPECIES)
            + ". A numeric NCBI taxonomy ID requires internet access."
        ),
    )
    optimization.add_argument(
        "--avoid-sequence",
        action="append",
        default=[],
        help="Additional DNA motif to avoid; may be passed multiple times.",
    )
    return parser


def config_from_args(args: argparse.Namespace) -> WorkflowConfig:
    default_avoid = ("AAAAA", "GGGGG", "CCCCC", "TTTTT")
    return WorkflowConfig(
        input_path=args.input,
        input_kind=args.input_kind,
        output_dir=args.output_dir,
        run_name=args.run_name,
        overhangs_path=args.overhangs,
        primers_path=args.primers,
        opool_length=args.opool_length,
        genes_per_subpool=args.genes_per_subpool,
        vector_oh1=args.vector_oh1,
        vector_oh2=args.vector_oh2,
        short_pool_max_size=args.short_pool_max_size,
        max_fragments=args.max_fragments,
        strip_nterm_met=not args.keep_nterm_met,
        typeiis_site=args.typeiis_site,
        typeiis_n=args.typeiis_n,
        additional_avoid_sequences=(*default_avoid, *args.avoid_sequence),
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        gc_window=args.gc_window,
        codon_species=args.codon_species,
        primer_mode=args.primer_mode,
        primer_start_at=args.primer_start_at,
        add_stuffer=not args.no_stuffer,
        random_seed=args.seed,
        skip_primer_assignment=args.skip_primers,
        exclude_unsafe_primers=not args.allow_unsafe_primers,
        show_progress=not args.quiet,
        force=args.force,
    )


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = run_workflow(config_from_args(args))
    except (ValueError, FileNotFoundError, FileExistsError, RuntimeError, AssertionError) as error:
        parser.exit(1, f"ERROR: {error}\n")
    return 0 if result.unassigned_genes == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
