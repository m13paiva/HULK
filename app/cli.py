#!/usr/bin/env python3
from pathlib import Path
import sys
import pandas as pd
import click

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

# -------- Pretty options: aligned columns + blank line between options --------
class SpacedCommand(click.Command):
    def format_options(self, ctx, formatter):
        rows = [p.get_help_record(ctx) for p in self.get_params(ctx)]
        rows = [r for r in rows if r]
        if not rows:
            return
        left_len = max(len(left) for left, _ in rows)
        pad = max(28, min(44, left_len + 2))
        width = formatter.width or 120
        with formatter.section("Options"):
            for i, (left, right) in enumerate(rows):
                right = right or ""
                formatter.write_text(f"  {left.ljust(pad)}{right}")
                if i < len(rows) - 1:
                    formatter.write_text("")

class HulkCommand(SpacedCommand):
    """Inject raw LOGO and raw HELP_BODY, preserving newlines exactly."""
    def get_help(self, ctx):
        base = super().get_help(ctx)

        # 1) Inject the raw ASCII LOGO right after "Usage:" (once)
        lines = base.splitlines()
        if lines and lines[0].startswith("Usage:"):
            rest = "\n".join(lines[1:])
            base = f"{lines[0]}\n\n{LOGO.rstrip()}\n\n{rest}"
        else:
            base = f"{LOGO.rstrip()}\n\n{base}"

        # 2) Insert raw HELP_BODY (verbatim) before "Options:"
        marker = "\nOptions:\n"
        if marker in base:
            head, tail = base.split(marker, 1)
            return f"{head}\n{HELP_BODY.rstrip()}\n\n{marker}{tail}"
        else:  # fallback
            return f"{base}\n\n{HELP_BODY.rstrip()}\n"

# Keep body text OUT of click's help (we inject it ourselves above)
HELP_BODY = (
    "HULK is a H(igh-volume b)ulk RNA-seq data preprocessing pipeline for NCBI SRA accessions.\n"
    "\n"
    "Given an input table (CSV/TSV/TXT) with columns 'Run', 'BioProject', and 'Model', HULK will:\n"
    "  • fetch SRA data (prefetch) with safe resume/overwrite,\n"
    "  • convert to FASTQ (fasterq-dump),\n"
    "  • perform QC/trimming with fastp,\n"
    "  • quantify against a transcriptome using kallisto.\n"
    "\n"
    "The transcriptome may be a compressed FASTA (.fa/.fasta/.fa.gz) or a prebuilt kallisto index (.idx).\n"
    "If a FASTA is provided, HULK builds a shared index once and reuses it.\n"
    "\n"
    "Concurrency is automatic: available CPU cores are detected (Docker/CGroups-aware), and work is split\n"
    "across SRAs with configurable per-job threading.\n"
    "\n"
    "Outputs are organized as <OUTPUT>/<BioProject>/<Run>/ ... A master log is written to <OUTPUT>/shared/log.txt.\n"
    "Re-runs are safe: completed SRAs are skipped unless --force is given; partial downloads are resumed/cleaned.\n"
    "\n"
    "Examples:\n"
    "  hulk -i samples.tsv -r transcripts.fasta.gz -o results\n"
    "  hulk -i samples.csv -r transcripts.idx --min-threads 2 --max-threads 8 -v\n"
    "  hulk -i table.txt -r species.fa.gz -t 8 -f -a\n"
)

@click.command(
    cls=HulkCommand,
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 180},
    help="",  # leave empty; we inject HELP_BODY verbatim ourselves
)
@click.version_option(__version__, "-V", "--version", prog_name="hulk")
@click.option("--smash", is_flag=True, hidden=True, is_eager=True, expose_value=False,
              callback=lambda ctx, p, v: (smash(), ctx.exit()) if v and not ctx.resilient_parsing else None)
# Required I/O
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Input table (.csv, .tsv, or .txt) with columns: 'Run', 'BioProject', 'Model'.",
)
@click.option(
    "-r", "--reference", "reference_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx).",
)
# Optional outputs
@click.option(
    "-o", "--output", "output_dir",
    type=click.Path(dir_okay=True, file_okay=False, path_type=Path),
    default=DEFAULT_OUTDIR, show_default=True,
    help="Output directory.",
)
# Performance
@click.option("--min-threads", type=int, default=DEFAULT_MIN_THREADS, show_default=True,
              help="Minimum number of threads per SRA.")
@click.option("-t", "--max-threads", type=int, default=DEFAULT_MAX_THREADS, show_default=True,
              help="Maximum total threads.")
# Flags
@click.option("-v", "--verbose", is_flag=True, help="Show live progress bars and console messages.")
@click.option("-f", "--force", "--overwrite", is_flag=True,
              help="Force re-run: overwrite totally/partially processed SRAs.")
@click.option("-a", "--aggregate", "--overall-table", is_flag=True,
              help="Create a merged TPM table across all BioProjects; if --gene-counts is set, also write a global gene-counts table.")
@click.option("-n", "--dry-run", is_flag=True,
              help="Validate inputs and configuration, print plan, and exit without running tools.")
@click.option(
    "-g", "--gene-counts", "tx2gene_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Enable gene counts using a tx2gene (.csv) with columns 'transcript_id','gene_id'.",
)
@click.option("--keep-fastq", is_flag=True, help="Keep FASTQ files.")
def cli(
    input_path: Path,
    reference_path: Path,
    output_dir: Path,
    min_threads: int,
    max_threads: int,
    verbose: bool,
    force: bool,
    aggregate: bool,
    dry_run: bool,
    tx2gene_path: Path | None,
    keep_fastq: bool,
):
    # Validate input extension
    suf = input_path.suffix.lower()
    if suf not in INPUT_FORMATS:
        raise click.UsageError(f"Input file must be one of: {', '.join(sorted(INPUT_FORMATS))} (got: {input_path.name})")

    # Read table
    try:
        if suf == ".csv":
            df = pd.read_csv(input_path, low_memory=False)
        elif suf == ".tsv":
            df = pd.read_csv(input_path, sep="\t", low_memory=False)
        else:
            df = pd.read_table(input_path)
    except Exception as e:
        raise click.ClickException(f"Could not read input file '{input_path}': {e}")

    if not REQUIRED_INPUT_COLS.issubset(df.columns):
        missing = REQUIRED_INPUT_COLS.difference(df.columns)
        raise click.UsageError(
            "Input file must contain columns: " + ", ".join(sorted(REQUIRED_INPUT_COLS)) +
            (f" (missing: {', '.join(sorted(missing))})" if missing else "")
        )
    df = df[list(REQUIRED_INPUT_COLS)]

    # Reference validation
    trans_suff = transcriptome_suffixes(reference_path)
    if trans_suff not in TRANSCRIPTOME_FORMATS:
        raise click.UsageError(
            f"Transcriptome file must be one of: {', '.join(sorted(TRANSCRIPTOME_FORMATS))} (got: {reference_path.name})"
        )

    # tx2gene validation (optional)
    if tx2gene_path is not None:
        try:
            tx2gene_df = pd.read_csv(tx2gene_path, sep=None, engine="python")
        except Exception as e:
            raise click.ClickException(f"Could not read tx2gene file '{tx2gene_path}': {e}")
        if not REQUIRED_TX2GENE_COLS.issubset(tx2gene_df.columns):
            raise click.UsageError(
                "tx2gene file must contain columns: " + ", ".join(sorted(REQUIRED_TX2GENE_COLS)) +
                f" (found: {', '.join(tx2gene_df.columns)})"
            )

    # Output dir
    try:
        output_dir = output_dir.expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise click.ClickException(f"Failed to prepare output directory '{output_dir}': {e}")

    # Dry-run
    if dry_run:
        click.secho("✅ Dry run: inputs and configuration validated.", fg="green")
        click.echo(f"- Input:        {input_path}")
        click.echo(f"- Reference:    {reference_path}")
        click.echo(f"- Output dir:   {output_dir}")
        click.echo(f"- Threads:      min={min_threads}, max={max_threads}")
        click.echo(f"- Force:        {force}")
        click.echo(f"- Aggregate:    {aggregate}")
        click.echo(f"- tx2gene:      {tx2gene_path if tx2gene_path else 'None'}")
        click.echo(f"- Keep FASTQ:   {keep_fastq}")
        sys.exit(0)

    # Run pipeline
    pipeline(
        df,
        output_dir,
        reference_path.expanduser().resolve(),
        min_threads,
        max_threads,
        verbose,
        force,
        keep_fastq,
        aggregate,
        tx2gene_path.expanduser().resolve() if tx2gene_path else None,
    )
    sys.exit(0)

def main():
    cli(standalone_mode=True)

if __name__ == "__main__":
    main()
