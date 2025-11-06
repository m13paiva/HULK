#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import json
from pathlib import Path
from typing import Any

import click
import pandas as pd
from click.core import ParameterSource

from .core import pipeline
from .entities import Config  # <-- use the Config class
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

CFG_ENV = "HULK_CONFIG"          # env var to point to a custom config path
CFG_DEFAULT_NAME = ".hulk.json"  # default config file in CWD
WIDE_HELP = 200                   # force very wide help for main + subcommands
RIGHT_COL = 50
LOGO_PAD = 20

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
    "  hulk -i samples.csv -r transcripts.idx --min-threads 2 --max-threads 8 -v\n"
    "  hulk -i table.txt -r species.fa.gz -t 8 -f -a\n"
)

# ---------------------------- Config helpers (persisted JSON for subcommands) ----------------------------
def _cfg_path(explicit: Path | None = None) -> Path:
    if explicit is not None:
        return explicit
    p_env = os.environ.get(CFG_ENV)
    if p_env:
        return Path(p_env).expanduser().resolve()
    # default to local file (works inside containers w/out $HOME perms)
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
    for k, v in payload.items():
        if v is not None:
            sect[k] = v
    cfg[section] = sect
    return _cfg_save(cfg, explicit)

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

            # fixed gap with minimum 2 spaces
            visible_left_len = len(left)
            gap = (RIGHT_COL - visible_left_len)
            if gap < 2:
                gap = 2

            formatter.write(f"  {left}{' ' * gap}{right}\n")

            if i < len(rows) - 1:
                formatter.write("\n")

# ---------------------------- Custom Group & Command ----------------------------
class HulkGroup(SpacedFormatterMixin, click.Group):
    """Custom Group that prints LOGO + HELP_BODY verbatim and formats sections."""
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
    """Subcommands use wide help + our spaced option formatting."""
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
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
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
              help="Minimum number of threads per SRR.")
@click.option("-t", "--max-threads", type=int, default=DEFAULT_MAX_THREADS, show_default=True,
              help="Maximum total threads.")
# Flags
@click.option("-v", "--verbose", is_flag=True, help="Show live progress bars and console messages.")
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
    verbose: bool,
    force: bool,
    aggregate: bool,
    dry_run: bool,
    tx2gene_path: Path | None,
    keep_fastq: bool,
):
    """
    Root group: run the pipeline when called directly; subcommands store defaults.
    """
    # If a subcommand (other than "run") is invoked, don't run pipeline.
    if ctx.invoked_subcommand is not None:
        return

    _run_pipeline(
        input_path=input_path,
        reference_path=reference_path,
        output_dir=output_dir,
        min_threads=min_threads,
        max_threads=max_threads,
        verbose=verbose,
        force=force,
        aggregate=aggregate,
        dry_run=dry_run,
        tx2gene_path=tx2gene_path,
        keep_fastq=keep_fastq,
    )

# --- Factor the actual execution into a helper used by both entry paths ---
def _run_pipeline(
    *,
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
) -> None:
    # Require inputs when running pipeline directly.
    if input_path is None or reference_path is None:
        raise click.UsageError("Missing required options: -i/--input and -r/--reference. See 'hulk -h'.")

    # Validate input table
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

    # Reference validation
    trans_suff = transcriptome_suffixes(reference_path)
    if trans_suff not in TRANSCRIPTOME_FORMATS:
        raise click.UsageError(
            f"Transcriptome file must be one of: {', '.join(sorted(TRANSCRIPTOME_FORMATS))} (got: {reference_path.name})"
        )

    # tx2gene validation
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

    # Build a Config object (this also creates <outdir>/shared and log.txt)
    cfg = Config(
        input_path=input_path,
        reference_path=reference_path,
        output_dir=output_dir,
        min_threads=min_threads,
        max_threads=max_threads,
        verbose=verbose,
        force=force,
        aggregate=aggregate,
        dry_run=dry_run,
        tx2gene_path=tx2gene_path,
        keep_fastq=keep_fastq,
    )

    # Dry run
    if dry_run:
        click.secho("\n✅ Dry run: inputs and configuration validated.\n", fg="green")
        click.echo(f"- Input:        {cfg.input_path}\n")
        click.echo(f"- Reference:    {cfg.reference_path}\n")
        click.echo(f"- Output dir:   {cfg.output_dir}\n")
        click.echo(f"- Threads:      min={cfg.min_threads}, max={cfg.max_threads}\n")
        click.echo(f"- Force:        {cfg.force}\n")
        click.echo(f"- Aggregate:    {cfg.aggregate}\n")
        click.echo(f"- tx2gene:      {cfg.tx2gene_path if cfg.tx2gene_path else 'None'}\n")
        click.echo(f"- Keep FASTQ:   {cfg.keep_fastq}\n")
        if cfg.trim_opts:
            click.echo(f"- Trim opts:    {cfg.trim_opts}")
        if cfg.tximport_opts:
            click.echo(f"- Tximport:     {cfg.tximport_opts}")
        sys.exit(0)

    # Run pipeline (keep current signature; supply params from cfg)
    pipeline(
        df,
        cfg.output_dir,
        cfg.reference_path.expanduser().resolve(),
        cfg.min_threads,
        cfg.max_threads,
        cfg.verbose,
        cfg.force,
        cfg.keep_fastq,
        cfg.aggregate,
        cfg.tx2gene_path.expanduser().resolve() if cfg.tx2gene_path else None,
        trim_opts=cfg.trim_opts or None,
        tximport_opts=cfg.tximport_opts or None,
    )
    sys.exit(0)

# ---------------------------- Subcommands ----------------------------
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
    help="Input table (.csv, .tsv, or .txt) with columns: 'Run', 'BioProject', 'Model'."
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
@click.option("-v", "--verbose", is_flag=True, help="Show live progress bars and console messages.")
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
    verbose: bool,
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
        verbose=verbose,
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
    "--ignore-tx-version",
    "ignore_tx_version",
    is_flag=True,
    default=False,
    help="If set, strip transcript version suffixes before matching (default: off).",
)
@click.option(
    "--config", "config_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=f"Path to settings JSON (default: ${CFG_ENV} or ./{CFG_DEFAULT_NAME})."
)
def tximport(mode: str | None, ignore_tx_version: bool, config_path: Path | None):
    """Save tximport defaults (used automatically by `hulk ...`)."""
    ctx = click.get_current_context()
    flag_provided = (ctx.get_parameter_source("ignore_tx_version") == ParameterSource.COMMANDLINE)

    if mode is None and not flag_provided:
        click.echo("Try: hulk tximport -h   (to see options)")
        return

    payload = {}
    if mode is not None:
        payload["mode"] = mode
    if flag_provided:
        payload["ignore_tx_version"] = ignore_tx_version  # True if flag used

    p = _cfg_update("tximport", payload, explicit=config_path)
    click.secho(f"Saved tximport settings to {p}", fg="green")

# ---------------------------- Entry point ----------------------------
def main():
    # Ensure wide formatter even if Click detects a narrow terminal
    os.environ.setdefault("COLUMNS", str(WIDE_HELP))
    cli(standalone_mode=True)

if __name__ == "__main__":
    main()
