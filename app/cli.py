#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import json
from pathlib import Path
from typing import Any, List

import click
import pandas as pd
from click.core import ParameterSource

from .core import pipeline                 # expects: pipeline(dataset, cfg)
from .entities import Config, Dataset
from .utils import transcriptome_suffixes, smash
from .logo import LOGO
from . import __version__

# ---------------------------- Constants ----------------------------
DEFAULT_OUTDIR = Path("output")
INPUT_FORMATS = {".csv", ".tsv", ".txt"}
REQUIRED_INPUT_COLS = {"Run", "BioProject", "Model"}
TRANSCRIPTOME_FORMATS = {".fasta.gz", ".fa.gz", ".fasta", ".fa", ".idx"}
DEFAULT_MIN_THREADS = 4
DEFAULT_MAX_THREADS = 10
REQUIRED_TX2GENE_COLS = {"transcript_id", "gene_id"}
CFG_FIXED_PATH = Path("/app/.hulk.json")
WIDE_HELP = 200                   # force very wide help for main + subcommands
RIGHT_COL = 50
LOGO_PAD = 20
HEAD_N = 5                        # rows/files to preview in "head"

# ---------------------------- Pretty help text (verbatim) ----------------------------
HELP_BODY = (
    "HULK is a H(igh-volume b)ulk RNA-seq data preprocessing pipeline for NCBI SRA accessions.\n"
    "\n"
    "Given an input table (CSV/TSV/TXT) with columns 'Run', 'BioProject', and 'Model', HULK will:\n"
    "  • fetch SRR data (prefetch) with safe resume/overwrite,\n"
    "  • convert to FASTQ (fasterq-dump),\n"
    "  • perform QC/trimming with fastp,\n"
    "  • quantify against a transcriptome using kallisto.\n"
    "\n"
    "The transcriptome may be a compressed FASTA (.fa/.fasta/.fa.gz) or a prebuilt kallisto index (.idx).\n"
    "If a FASTA is provided, HULK builds a shared index once and reuses it.\n"
    "\n"
    "Concurrency is automatic: available CPU cores are detected (Docker/CGroups-aware), and work is split\n"
    "across SRRs with configurable per-job threading.\n"
    "\n"
    "Outputs are organized as <OUTPUT>/<BioProject>/<Run>/ ... A master log is written to <OUTPUT>/shared/log.txt.\n"
    "Re-runs are safe: completed SRRs are skipped unless --force is given; partial downloads are resumed/cleaned.\n"
    "\n"
    "Examples:\n"
    "  hulk -i samples.tsv -r transcripts.fasta.gz -o results\n"
    "  hulk -i samples.csv -r transcripts.idx --min-threads 2 --max-threads 8 --no-verbosity\n"
    "  hulk -i reads_dir/ -r species.fa.gz -t 8 -f -a -y\n"
)

# ---------------------------- Config helpers (persisted JSON for subcommands) ----------------------------

def _cfg_path(_: Path | None = None) -> Path:
    # Always use the fixed path
    return CFG_FIXED_PATH

def _cfg_load(_: Path | None = None) -> dict[str, Any]:
    p = _cfg_path(None)
    if not p.exists():
        return {}
    try:
        with open(p, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}

def _cfg_save(cfg: dict[str, Any], _: Path | None = None) -> Path:
    p = _cfg_path(None)
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    tmp = p.with_suffix(p.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    tmp.replace(p)
    return p

def _cfg_update(section: str, payload: dict[str, Any], _: Path | None = None) -> Path:
    cfg = _cfg_load(None)
    sect = cfg.get(section, {})
    for k, v in payload.items():
        if v is not None:
            sect[k] = v
    cfg[section] = sect
    return _cfg_save(cfg, None)

def _cfg_reset(section: str, _: Path | None = None) -> Path:
    cfg = _cfg_load(None)
    if section in cfg:
        del cfg[section]
    return _cfg_save(cfg, None)
# ---------------------------- Pretty formatting mixin ----------------------------
def _pad_block(text: str, n: int) -> str:
    lines = text.rstrip("\n").splitlines()
    return "\n".join((" " * n) + line for line in lines)

class SpacedFormatterMixin:
    """Align left/right columns and insert blank lines between rows, preserving text verbatim."""
    def _fmt_rows_with_spacing(self, formatter: click.HelpFormatter, title: str, rows: list[tuple[str, str | None]]) -> None:
        rows = [r for r in rows if r]
        if not rows:
            return
        formatter.write(f"{title}:\n")
        for i, (left, right) in enumerate(rows):
            left = left or ""
            right = right or ""
            visible_left_len = len(left)
            gap = (RIGHT_COL - visible_left_len)
            if gap < 2:
                gap = 2
            formatter.write(f"  {left}{' ' * gap}{right}\n")
            if i < len(rows) - 1:
                formatter.write("\n")

# ---------------------------- Custom Group & Command ----------------------------
class HulkGroup(SpacedFormatterMixin, click.Group):
    def make_context(self, info_name, args, parent=None, **extra):
        ctx = super().make_context(info_name, args, parent=parent, **extra)
        ctx.max_content_width = WIDE_HELP
        return ctx

    def format_help(self, ctx, formatter):
        self.format_usage(ctx, formatter)
        formatter.write_text("\n")
        formatter.write(_pad_block(LOGO, LOGO_PAD) + "\n\n")
        formatter.write(HELP_BODY.rstrip() + "\n\n")
        self.format_commands(ctx, formatter)
        formatter.write("\n")
        self.format_options(ctx, formatter)

    def format_options(self, ctx, formatter):
        opts = [p.get_help_record(ctx) for p in self.get_params(ctx)]
        opts = [o for o in opts if o]
        self._fmt_rows_with_spacing(formatter, "Options", opts)

    def format_commands(self, ctx, formatter):
        commands: list[tuple[str, str]] = []
        for name in self.list_commands(ctx):
            cmd = self.get_command(ctx, name)
            if not cmd or cmd.hidden:
                continue
            one_liner = (cmd.help or "").strip().splitlines()[0]
            commands.append((name, one_liner))
        if commands:
            self._fmt_rows_with_spacing(formatter, "Commands", commands)

class HulkCommand(SpacedFormatterMixin, click.Command):
    def make_context(self, info_name, args, parent=None, **extra):
        ctx = super().make_context(info_name, args, parent=parent, **extra)
        ctx.max_content_width = WIDE_HELP
        return ctx

    def format_help(self, ctx, formatter):
        self.format_usage(ctx, formatter)
        formatter.write("\n")
        if self.help:
            formatter.write(self.help.strip() + "\n\n")
        self.format_options(ctx, formatter)

    def format_options(self, ctx, formatter):
        opts = [p.get_help_record(ctx) for p in self.get_params(ctx)]
        opts = [o for o in opts if o]
        self._fmt_rows_with_spacing(formatter, "Options", opts)

# ---------------------------- Root CLI (Group) ----------------------------
@click.group(
    cls=HulkGroup,
    invoke_without_command=True,
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": WIDE_HELP},
)
@click.version_option(__version__, "-V", "--version", prog_name="hulk")
@click.option("--smash", is_flag=True, hidden=True, is_eager=True, expose_value=False,
              callback=lambda ctx, p, v: (smash(), ctx.exit()) if v and not ctx.resilient_parsing else None)
# Required I/O (optional here so subcommands can run alone)
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),  # dir_okay=True for FASTQ mode
    required=False,
    help="Input table (.csv/.tsv/.txt) with 'Run','BioProject','Model' OR a directory of FASTQ files.",
)
@click.option(
    "-r", "--reference", "reference_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
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
              help="Minimum number of threads per SRR.")
@click.option("-t", "--max-threads", type=int, default=DEFAULT_MAX_THREADS, show_default=True,
              help="Maximum total threads.")
# Flags
@click.option("--verbosity/--no-verbosity", default=True, show_default=True,
              help="Show live progress bars and console messages (default: on).")
@click.option("-y", "--yes", is_flag=True, help="Assume 'yes' to prompts and run without asking.")
@click.option("-f", "--force", "--overwrite", is_flag=True,
              help="Force re-run: overwrite totally/partially processed SRRs.")
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
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path | None,
    reference_path: Path | None,
    output_dir: Path,
    min_threads: int,
    max_threads: int,
    verbosity: bool,
    yes: bool,
    force: bool,
    aggregate: bool,
    dry_run: bool,
    tx2gene_path: Path | None,
    keep_fastq: bool,
):
    """Root group: run the pipeline when called directly; subcommands store defaults."""
    if ctx.invoked_subcommand is not None:
        return
    _run_pipeline(
        input_path=input_path,
        reference_path=reference_path,
        output_dir=output_dir,
        min_threads=min_threads,
        max_threads=max_threads,
        verbosity=verbosity,
        yes=yes,
        force=force,
        aggregate=aggregate,
        dry_run=dry_run,
        tx2gene_path=tx2gene_path,
        keep_fastq=keep_fastq,
    )

# ---------------------------- Helpers ----------------------------
def _print_plan_for_srr(df: pd.DataFrame, dataset: Dataset, cfg: Config) -> None:
    bps = sorted(dataset.bioprojects.keys())
    n_bps = len(bps)
    n_samples = len(dataset.samples)
    click.echo(f"\nDetected SRR run:")
    click.echo(f"- #SRR samples: {n_samples}")
    click.echo(f"- #BioProjects: {n_bps}")
    if n_bps:
        click.echo(f"- BioProjects:  {', '.join(bps)}")
    click.echo("\nInput table (head):")
    click.echo(df.head(HEAD_N).to_string(index=False))
    click.echo("\nConfig:")
    click.echo(cfg.summary())

def _scan_fastqs(directory: Path) -> List[Path]:
    exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    files = sorted(p for p in Path(directory).iterdir() if p.is_file() and p.name.lower().endswith(exts))
    return files

def _print_plan_for_fastq(fastq_files: List[Path], dataset: Dataset, cfg: Config) -> None:
    click.echo(f"\nDetected FASTQ run:")
    click.echo(f"- #FASTQ files (samples): {len(fastq_files)}")
    click.echo(f"- #BioProjects: N/A")
    head = [f.name for f in fastq_files[:HEAD_N]]
    click.echo("\nFASTQ files (head):")
    for name in head:
        click.echo(f"  - {name}")
    click.echo("\nConfig:")
    click.echo(cfg.summary())

# ---------------------------- Runner ----------------------------
def _run_pipeline(
    *,
    input_path: Path | None,
    reference_path: Path | None,
    output_dir: Path,
    min_threads: int,
    max_threads: int,
    verbosity: bool,
    yes: bool,
    force: bool,
    aggregate: bool,
    dry_run: bool,
    tx2gene_path: Path | None,
    keep_fastq: bool,
) -> None:
    # Require inputs
    if input_path is None or reference_path is None:
        raise click.UsageError("Missing required options: -i/--input and -r/--reference. See 'hulk -h'.")

    # Build Config (creates <outdir>/shared and log.txt). NOTE: Config now expects 'outdir'.
    cfg = Config(
        input_path=input_path,
        reference_path=reference_path,
        outdir=output_dir,                  # <-- changed name
        min_threads=min_threads,
        max_threads=max_threads,
        verbose=verbosity,
        force=force,
        aggregate=aggregate,
        dry_run=dry_run,
        tx2gene=tx2gene_path,
        keep_fastq=keep_fastq,
    )

    # Validate reference format (common to both modes)
    trans_suff = transcriptome_suffixes(reference_path)
    if trans_suff not in TRANSCRIPTOME_FORMATS:
        raise click.UsageError(
            f"Transcriptome file must be one of: {', '.join(sorted(TRANSCRIPTOME_FORMATS))} (got: {reference_path.name})"
        )

    # Build Dataset depending on input_path type
    if input_path.is_dir():
        # FASTQ mode
        fastq_files = _scan_fastqs(input_path)
        if not fastq_files:
            raise click.UsageError(f"No FASTQ files found in directory: {input_path}")
        dataset = Dataset.from_fastq_dir(input_path, cfg)
        _print_plan_for_fastq(fastq_files, dataset, cfg)
        df = None
    else:
        # SRR mode: validate table then build Dataset
        suf = input_path.suffix.lower()
        if suf not in INPUT_FORMATS:
            raise click.UsageError(f"Input file must be one of: {', '.join(sorted(INPUT_FORMATS))} (got: {input_path.name})")
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
        dataset = Dataset.from_dataframe(df, cfg)
        _print_plan_for_srr(df, dataset, cfg)

    # tx2gene validation (only if provided)
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

    # If dry run, stop after printing the plan
    if dry_run:
        click.secho("\n✅ Dry run complete. No tools executed.\n", fg="green")
        sys.exit(0)

    # Ask for confirmation unless -y/--yes was provided
    if not yes:
        click.echo("")  # spacing
        click.confirm("Proceed with the run? (Y/n)", default=True, abort=True)

    # Run the pipeline
    pipeline(dataset, cfg)
    sys.exit(0)

# ---------------------------- Subcommands ----------------------------
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    required=False,
    help="Input table (.csv/.tsv/.txt) with 'Run','BioProject','Model' OR a directory of FASTQ files."
)
@click.option(
    "-r", "--reference", "reference_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
    help="Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx)."
)
@click.option(
    "-o", "--output", "output_dir",
    type=click.Path(dir_okay=True, file_okay=False, path_type=Path),
    default=DEFAULT_OUTDIR, show_default=True,
    help="Output directory."
)
@click.option("--min-threads", type=int, default=DEFAULT_MIN_THREADS, show_default=True,
              help="Minimum number of threads per SRR.")
@click.option("-t", "--max-threads", type=int, default=DEFAULT_MAX_THREADS, show_default=True,
              help="Maximum total threads.")
@click.option("--verbosity/--no-verbosity", default=True, show_default=True,
              help="Show live progress bars and console messages (default: on).")
@click.option("-y", "--yes", is_flag=True, help="Assume 'yes' to prompts and run without asking.")
@click.option("-f", "--force", "--overwrite", is_flag=True,
              help="Force re-run: overwrite totally/partially processed SRRs.")
@click.option("-a", "--aggregate", "--overall-table", is_flag=True,
              help="Create a merged TPM table; with --gene-counts, also a global gene-counts table.")
@click.option("-n", "--dry-run", is_flag=True,
              help="Validate inputs and configuration, print plan, and exit without running tools.")
@click.option(
    "-g", "--gene-counts", "tx2gene_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Enable gene counts using a tx2gene (.csv) with columns 'transcript_id','gene_id'."
)
@click.option("--keep-fastq", is_flag=True, help="Keep FASTQ files.")
def run_cmd(
    input_path: Path | None,
    reference_path: Path | None,
    output_dir: Path,
    min_threads: int,
    max_threads: int,
    verbosity: bool,
    yes: bool,
    force: bool,
    aggregate: bool,
    dry_run: bool,
    tx2gene_path: Path | None,
    keep_fastq: bool,
):
    """Entry-point that mirrors the root options."""
    _run_pipeline(
        input_path=input_path,
        reference_path=reference_path,
        output_dir=output_dir,
        min_threads=min_threads,
        max_threads=max_threads,
        verbosity=verbosity,
        yes=yes,
        force=force,
        aggregate=aggregate,
        dry_run=dry_run,
        tx2gene_path=tx2gene_path,
        keep_fastq=keep_fastq,
    )

@cli.command("trim", cls=HulkCommand, help="Configure fastp trimming defaults.")
@click.option("-ws", "--window-size", type=int, default=None,
              help="fastp sliding window size (integer).")
@click.option("-mq", "--mean-quality", type=int, default=None,
              help="fastp mean quality threshold (integer).")
@click.option("--reset-defaults", is_flag=True, help="Reset trim options to built-in defaults.")
def trim(window_size: int | None, mean_quality: int | None, reset_defaults: bool):
    """Save trimming defaults (used automatically by `hulk ...`)."""
    if reset_defaults:
        p = _cfg_reset("trim")
        click.secho(f"Reset trim settings to defaults at {p}", fg="green")
        return

    if window_size is None and mean_quality is None:
        click.echo("Try: hulk trim -h   (to see options)")
        return

    p = _cfg_update("trim", {"window_size": window_size, "mean_quality": mean_quality})
    click.secho(f"Saved trim settings to {p}", fg="green")

@cli.command("tximport", cls=HulkCommand, help="Configure tximport aggregation/normalization.")
@click.option(
    "-m", "--mode",
    type=click.Choice(
        ["raw_counts", "length_scaled_tpm", "scaled_tpm", "dtu_scaled_tpm"],
        case_sensitive=False
    ),
    default=None,
    help="tximport aggregation/normalization mode.",
)
@click.option(
    "--ignore-tx-version",
    "ignore_tx_version",
    is_flag=True,
    default=False,
    help="If set, strip transcript version suffixes before matching (default: off).",
)
@click.option("--reset-defaults", is_flag=True, help="Reset tximport options to built-in defaults.")
def tximport(mode: str | None, ignore_tx_version: bool, reset_defaults: bool):
    """Save tximport defaults (used automatically by `hulk ...`)."""
    if reset_defaults:
        p = _cfg_reset("tximport")
        click.secho(f"Reset tximport settings to defaults at {p}", fg="green")
        return

    # Only persist the flag if the user actually passed it on the command line
    ctx = click.get_current_context()
    flag_provided = (ctx.get_parameter_source("ignore_tx_version") == ParameterSource.COMMANDLINE)

    if mode is None and not flag_provided:
        click.echo("Try: hulk tximport -h   (to see options)")
        return

    payload = {}
    if mode is not None:
        payload["mode"] = mode
    if flag_provided:
        payload["ignore_tx_version"] = ignore_tx_version

    p = _cfg_update("tximport", payload)
    click.secho(f"Saved tximport settings to {p}", fg="green")

# ---------------------------- Entry point ----------------------------
def main():
    os.environ.setdefault("COLUMNS", str(WIDE_HELP))
    cli(standalone_mode=True)

if __name__ == "__main__":
    main()
