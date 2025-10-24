#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import json
from pathlib import Path
from typing import Any

import click
import pandas as pd

from .core import pipeline
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

CFG_ENV = "HULK_CONFIG"  # env var to point to a custom config path
CFG_DEFAULT_NAME = ".hulk.json"

# ---------------------------- Pretty help ----------------------------
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

class SpacedFormatterMixin:
    """Mix-in to format options/commands with aligned columns + blank lines."""
    def _fmt_rows_with_spacing(self, formatter, title: str, rows: list[tuple[str, str | None]]) -> None:
        rows = [r for r in rows if r]
        if not rows:
            return
        left_len = max(len(left) for left, _ in rows)
        pad = max(28, min(44, left_len + 2))
        with formatter.section(title):
            for i, (left, right) in enumerate(rows):
                formatter.write_text(f"  {left.ljust(pad)}{(right or '')}")
                if i < len(rows) - 1:
                    formatter.write_text("")

class HulkGroup(SpacedFormatterMixin, click.Group):
    """Custom Group to inject LOGO, body text, and pretty sections."""
    def format_help(self, ctx, formatter):
        # Print "Usage:" line first (Click default)
        self.format_usage(ctx, formatter)
        formatter.write_text("")  # blank line

        # Inject raw LOGO exactly once
        formatter.write_text(LOGO.rstrip())
        formatter.write_text("")  # blank line

        # Inject the body paragraphs verbatim
        formatter.write_text(HELP_BODY.rstrip())
        formatter.write_text("")  # blank line

        # Then show commands and options in our pretty formatters
        self.format_commands(ctx, formatter)
        self.format_options(ctx, formatter)

    def format_options(self, ctx, formatter):
        # Collect normal options for this group (not subcommands)
        opts = [p.get_help_record(ctx) for p in self.get_params(ctx)]
        opts = [o for o in opts if o]
        self._fmt_rows_with_spacing(formatter, "Options", opts)

    def format_commands(self, ctx, formatter):
        # Show subcommands (+ a one-line description each)
        cmds: dict[str, click.Command] = self.list_commands(ctx)
        if not cmds:
            return
        rows: list[tuple[str, str]] = []
        for name in cmds:
            cmd = self.get_command(ctx, name)
            if cmd is None or cmd.hidden:
                continue
            help_line = (cmd.help or "").strip().splitlines()[0]
            rows.append((name, help_line))
        if rows:
            self._fmt_rows_with_spacing(formatter, "Commands", rows)

class HulkCommand(SpacedFormatterMixin, click.Command):
    """For subcommands to inherit our option spacing (logo/body only on group)."""
    def format_help(self, ctx, formatter):
        # Standard Click header for subcommands
        self.format_usage(ctx, formatter)
        formatter.write_text("")  # blank line
        # One-liner description if present
        if self.help:
            formatter.write_text(self.help.strip())
            formatter.write_text("")
        # Options in spaced style
        self.format_options(ctx, formatter)

# ---------------------------- Config helpers ----------------------------
def _cfg_path(explicit: Path | None = None) -> Path:
    if explicit is not None:
        return explicit
    p_env = os.environ.get(CFG_ENV)
    if p_env:
        return Path(p_env).expanduser().resolve()
    # default to a local file in CWD (works in containers without $HOME perms)
    return Path.cwd() / CFG_DEFAULT_NAME

def _cfg_load(explicit: Path | None = None) -> dict[str, Any]:
    p = _cfg_path(explicit)
    if not p.exists():
        return {}
    try:
        with open(p, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}

def _cfg_save(cfg: dict[str, Any], explicit: Path | None = None) -> Path:
    p = _cfg_path(explicit)
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        # If parent can't be created (e.g., CWD), we still try to write the file
        pass
    tmp = p.with_suffix(p.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    tmp.replace(p)
    return p

def _cfg_update(section: str, payload: dict[str, Any], explicit: Path | None = None) -> Path:
    cfg = _cfg_load(explicit)
    sect = cfg.get(section, {})
    # only set provided (non-None) keys
    for k, v in payload.items():
        if v is not None:
            sect[k] = v
    cfg[section] = sect
    return _cfg_save(cfg, explicit)

# ---------------------------- Root CLI (Group) ----------------------------
@click.group(
    cls=HulkGroup,
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 180},
)
@click.version_option(__version__, "-V", "--version", prog_name="hulk")
@click.option("--smash", is_flag=True, hidden=True, is_eager=True, expose_value=False,
              callback=lambda ctx, p, v: (smash(), ctx.exit()) if v and not ctx.resilient_parsing else None)
# Required I/O
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,  # allow calling subcommands without these
    help="Input table (.csv, .tsv, or .txt) with columns: 'Run', 'BioProject', 'Model'.",
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
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path | None,
    reference_path: Path | None,
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
    """
    Root group: if called with options (and not a subcommand), run the pipeline.
    Subcommands (e.g. `trim`, `tximport`) configure defaults saved to a JSON file.
    """
    # If a subcommand was invoked, just stash the common options and return.
    if ctx.invoked_subcommand is not None:
        # Stash for subcommand usage if needed
        ctx.obj = {
            "output_dir": output_dir,
        }
        return

    # If no subcommand: we expect the pipeline arguments to be present.
    if input_path is None or reference_path is None:
        raise click.UsageError("Missing required options: -i/--input and -r/--reference. See 'hulk -h'.")

    # Validate input table extension
    suf = input_path.suffix.lower()
    if suf not in INPUT_FORMATS:
        raise click.UsageError(
            f"Input file must be one of: {', '.join(sorted(INPUT_FORMATS))} (got: {input_path.name})"
        )

    # Read sample table
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

    # Load saved subcommand settings and merge into kwargs for pipeline
    cfg = _cfg_load()
    trim_opts = cfg.get("trim") or {}
    tximport_opts = cfg.get("tximport") or {}

    # Dry-run summary
    if dry_run:
        click.secho("\n✅ Dry run: inputs and configuration validated.\n", fg="green")
        click.echo(f"- Input:        {input_path}\n")
        click.echo(f"- Reference:    {reference_path}\n")
        click.echo(f"- Output dir:   {output_dir}\n")
        click.echo(f"- Threads:      min={min_threads}, max={max_threads}\n")
        click.echo(f"- Force:        {force}\n")
        click.echo(f"- Aggregate:    {aggregate}\n")
        click.echo(f"- tx2gene:      {tx2gene_path if tx2gene_path else 'None'}\n")
        click.echo(f"- Keep FASTQ:   {keep_fastq}\n")
        if trim_opts:
            click.echo(f"- Trim opts:    {trim_opts}")
        if tximport_opts:
            click.echo(f"- Tximport:     {tximport_opts}")
        sys.exit(0)

    # Run pipeline — ensure it accepts trim_opts / tximport_opts (add in core if needed)
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
        trim_opts=trim_opts or None,
        tximport_opts=tximport_opts or None,
    )
    sys.exit(0)

# ---------------------------- Subcommands ----------------------------
@cli.command("trim", cls=HulkCommand, help="Configure fastp trimming defaults.")
@click.option("-ws", "--window-size", type=int, default=None,
              help="fastp sliding window size (integer).")
@click.option("-mq", "--mean-quality", type=int, default=None,
              help="fastp mean quality threshold (integer).")
@click.option("--config", "config_path", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help=f"Path to settings JSON (default: ${CFG_ENV} or ./{CFG_DEFAULT_NAME}).")
def trim(window_size: int | None, mean_quality: int | None, config_path: Path | None):
    """Save trimming defaults (used automatically by `hulk ...`)."""
    if window_size is None and mean_quality is None:
        click.echo("Try: hulk trim -h   (to see options)")
        return
    p = _cfg_update("trim",
                    {"window_size": window_size, "mean_quality": mean_quality},
                    explicit=config_path)
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
    "--ignore_tx_version/--no-ignore_tx_version",
    default=None,
    help="Strip transcript version suffixes before matching (default: false).",
)
@click.option("--config", "config_path", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help=f"Path to settings JSON (default: ${CFG_ENV} or ./{CFG_DEFAULT_NAME}).")
def tximport(mode: str | None, ignore_tx_version: bool | None, config_path: Path | None):
    """Save tximport defaults (used automatically by `hulk ...`)."""
    if mode is None and ignore_tx_version is None:
        click.echo("Try: hulk tximport -h   (to see options)")
        return
    p = _cfg_update("tximport",
                    {"mode": mode, "ignore_tx_version": ignore_tx_version},
                    explicit=config_path)
    click.secho(f"Saved tximport settings to {p}", fg="green")

# ---------------------------- Entry point ----------------------------
def main():
    cli(standalone_mode=True)

if __name__ == "__main__":
    main()
