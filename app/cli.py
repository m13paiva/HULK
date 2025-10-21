#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import pandas as pd
from .core import pipeline
from .utils import transcriptome_suffixes, smash
from .logo import LOGO
from . import __version__

DEFAULT_OUTDIR = Path("output")
INPUT_FORMATS = {".csv", ".tsv", ".txt"}
REQUIRED_INPUT_COLS = {"Run", "BioProject", "Model"}
TRANSCRIPTOME_FORMATS = {".fasta.gz", ".fa.gz", ".fasta", ".fa", ".idx"}
DEFAULT_MIN_THREADS = 4
DEFAULT_MAX_THREADS = 10
REQUIRED_TX2GENE_COLS = {"transcript_id", "gene_id"}


def parse_args():
    # --- lightweight pre-parser ONLY for --smash (doesn't steal -h/--help)
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--smash", action="store_true")
    ns_pre, _ = pre.parse_known_args()
    if ns_pre.smash:
        smash()
        sys.exit(0)

    # --- full parser
    parser = argparse.ArgumentParser(
        prog="hulk",
        description=(
            LOGO
            + "HULK is a H(igh-volume b)ulk RNA-seq data preprocessing pipeline for NCBI SRA accessions.\n\n"
              "Given an input table (CSV/TSV/TXT) with columns 'Run', 'BioProject', and 'Model', "
              "HULK will:\n"
              "  • fetch SRA data (prefetch) with safe resume/overwrite,\n"
              "  • convert to FASTQ (fasterq-dump),\n"
              "  • perform QC/trimming with fastp,\n"
              "  • quantify against a transcriptome using kallisto.\n\n"
              "The transcriptome may be a compressed FASTA (.fa/.fasta/.fna.gz) or a \n"
              "prebuilt kallisto index (.idx). If a FASTA is provided, HULK builds a \n"
              "shared index once per run and reuses it.\n\n"
              "Concurrency is automatic: available CPU cores are detected (Docker/CGroups-aware), \n"
              "and work is split across SRAs with configurable per-job threading.\n"
        ),
        epilog=(
            "Outputs are organized as <OUTPUT>/<BioProject>/<Run>/ ...\n"
            "A master log is written to <OUTPUT>/shared/log.txt capturing all tool output.\n"
            "Re-runs are safe: completed SRAs are skipped unless --force is given; \n"
            "partial downloads are resumed/cleaned automatically.\n\n"
            "Examples:\n"
            "  hulk -i samples.tsv -r transcripts.fasta.gz -o results\n"
            "  hulk -i samples.csv -r transcripts.idx --min-threads 2 --max-threads 8 -v\n"
            "  hulk -i table.txt -r species.fa.gz -t 8 -f -a\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )

    # Version
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show version and exit."
    )

    # Hidden easter egg remains available, but not parsed before help anymore
    parser.add_argument("--smash", action="store_true", help=argparse.SUPPRESS)

    # Required I/O
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input table (.csv, .tsv, or .txt) with columns: 'Run', 'BioProject', 'Model'."
    )
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx)."
    )

    # Optional outputs
    parser.add_argument(
        "-o", "--output",
        default=Path(DEFAULT_OUTDIR),
        help="Output directory (default: 'output')."
    )

    # Performance
    parser.add_argument(
        "--min-threads",
        type=int,
        default=DEFAULT_MIN_THREADS,
        help=f"Minimum number of threads per SRA (default: {DEFAULT_MIN_THREADS})."
    )
    parser.add_argument(
        "-t", "--max-threads",
        type=int,
        default=DEFAULT_MAX_THREADS,
        help=f"Maximum total threads (default: {DEFAULT_MAX_THREADS})."
    )

    # Flags
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Show live progress bars and console messages."
    )

    parser.add_argument(
        "-f", "--force",
        action="store_true",
        help="Force re-run: overwrite totally/partially processed SRAs."
    )
    parser.add_argument(
        "--overwrite",
        dest="force",
        action="store_true",
        help=argparse.SUPPRESS
    )

    parser.add_argument(
        "-a", "--aggregate",
        action="store_true",
        help="Create a merged TPM table across all BioProjects; "
             "if --gene-counts is set, also write a global gene-counts table."
    )
    parser.add_argument(
        "--overall-table",
        dest="aggregate",
        action="store_true",
        help=argparse.SUPPRESS
    )
    """to do
    parser.add_argument(
        "-n", "--dry-run",
        action="store_true",
        help="Plan and validate without executing external tools."
    )
    """

    parser.add_argument(
        "-g", "--gene-counts",
        metavar="TX2GENE_FILE",
        default=None,
        help="Enable gene counts using a tx2gene (.csv) with columns 'transcript_id','gene_id'."
    )
    parser.add_argument(
        "--keep-fastq",
        action="store_true",
        help="Keep FASTQ files."
    )

    args = parser.parse_args()

    # --- validation: input file exists ---
    input_path = Path(args.input)
    if not input_path.is_file():
        parser.error(f"Input file does not exist: {input_path}")

    # --- validation: correct extension ---
    if input_path.suffix.lower() not in INPUT_FORMATS:
        parser.error(f"Input file must be one of: {', '.join(sorted(INPUT_FORMATS))}")

    # --- validation: required columns ---
    try:
        if input_path.suffix.lower() == ".csv":
            df = pd.read_csv(input_path, low_memory=False)
        elif input_path.suffix.lower() == ".tsv":
            df = pd.read_csv(input_path, sep="\t", low_memory=False)
        else:
            df = pd.read_table(input_path)
    except Exception as e:
        parser.error(f"Could not read input file: {e}")

    if not REQUIRED_INPUT_COLS.issubset(df.columns):
        parser.error(f"Input file must contain columns: {', '.join(sorted(REQUIRED_INPUT_COLS))}")

    df = df[list(REQUIRED_INPUT_COLS)]
    args.df = df

    # --- transcriptome / index validation ---
    reference_path = Path(args.reference).expanduser().resolve()
    if not reference_path.is_file():
        parser.error(f"reference file does not exist: {reference_path}")

    trans_suff = transcriptome_suffixes(reference_path)
    if trans_suff not in TRANSCRIPTOME_FORMATS:
        parser.error(f"Transcriptome file must be one of: {', '.join(sorted(TRANSCRIPTOME_FORMATS))}")

    args.reference = reference_path

    # --- tx2gene validation (optional) ---
    if args.gene_counts is not None:
        tx2gene_path = Path(args.gene_counts).expanduser().resolve()
        if not tx2gene_path.is_file():
            parser.error(f"tx2gene file does not exist: {tx2gene_path}")
        try:
            tx2gene_df = pd.read_csv(tx2gene_path, sep=None, engine="python")
        except Exception as e:
            parser.error(f"Could not read tx2gene file: {tx2gene_path} ({e})")
        if not REQUIRED_TX2GENE_COLS.issubset(tx2gene_df.columns):
            parser.error(
                "tx2gene file must contain columns: "
                + ", ".join(sorted(REQUIRED_TX2GENE_COLS))
                + f" (found: {', '.join(tx2gene_df.columns)})"
            )
        args.gene_counts = tx2gene_path

    # --- output dir normalization ---
    if args.output:
        args.output = Path(args.output).expanduser().resolve()
        if not args.output.is_dir():
            args.output.mkdir(parents=True)

    return args


def main():
    args = parse_args()
    pipeline(
        args.df,
        args.output,
        args.reference,
        args.min_threads,
        args.max_threads,
        args.verbose,
        args.force,
        args.keep_fastq,
        args.aggregate,
        args.gene_counts,
        # NOTE: args.dry_run is parsed; wire it into pipeline when implemented
    )


if __name__ == "__main__":
    main()
