#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import json
from pathlib import Path
from typing import Any, List, Dict, Optional

import click
import pandas as pd
from click.core import ParameterSource

from .core import pipeline                 # expects: pipeline(dataset, cfg)
from .entities import Config, Dataset
from .utils import transcriptome_suffixes, smash, scan_fastqs
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
WIDE_HELP = 200
RIGHT_COL = 50
LOGO_PAD = 20

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
    "Outputs are organized as <OUTPUT>/<BioProject>/<Run>/ ... A master log is written to <OUTPUT>/shared/log.txt.\n"
    "Re-runs are safe: completed SRRs are skipped unless --force is given; partial downloads are resumed/cleaned.\n"
)

# ---------------------------- Persisted config I/O ----------------------------

def _cfg_path(_: Path | None = None) -> Path:
    return CFG_FIXED_PATH

def _cfg_load(_: Path | None = None) -> Dict[str, Any]:
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
    p.parent.mkdir(parents=True, exist_ok=True)
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

# ---------------------------- Help formatting ----------------------------

def _pad_block(text: str, n: int) -> str:
    lines = text.rstrip("\n").splitlines()
    return "\n".join((" " * n) + line for line in lines)

class SpacedFormatterMixin:
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

# ---------------------------- Root CLI ----------------------------

@click.group(
    cls=HulkGroup,
    invoke_without_command=True,
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": WIDE_HELP},
)
@click.version_option(__version__, "-V", "--version", prog_name="hulk")
@click.option("--smash", is_flag=True, hidden=True, is_eager=True, expose_value=False,
              callback=lambda ctx, p, v: (smash(), ctx.exit()) if v and not ctx.resilient_parsing else None)
@click.option("-i", "--input", "input_path",
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=False,
              help="Input table (.csv/.tsv/.txt) with 'Run','BioProject','Model' OR a directory of FASTQ files.")
@click.option("-r", "--reference", "reference_path",
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              required=False,
              help="Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx).")
@click.option("-o", "--output", "output_dir",
              type=click.Path(dir_okay=True, file_okay=False, path_type=Path),
              default=DEFAULT_OUTDIR, show_default=True,
              help="Output directory.")
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
              help="Create a merged TPM table across all BioProjects; with --gene-counts, also a global gene-counts table.")
@click.option("-n", "--dry-run", is_flag=True,
              help="Validate inputs and configuration, print plan, and exit without running tools.")
@click.option("-g", "--gene-counts", "tx2gene_path",
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None,
              help="Enable gene counts using a tx2gene (.csv) with columns 'transcript_id','gene_id'.")
@click.option(
    "--no-cache",
    is_flag=True,
    help="Disable SRA caching (do not allocate a shared cache volume).",
)
@click.option(
    "-c", "--cache",
    "cache_gb",
    type=int,
    default=None,
    help="Maximum SRA cache size (GiB). Default: auto (300 GiB or free space based).",
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
    no_cache: bool,
    cache_gb: int | None,
    keep_fastq: bool,
):
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
        no_cache=no_cache,
        cache_gb=cache_gb,
        keep_fastq=keep_fastq,
    )

# ---------------------------- Summary printing ----------------------------

def _val(x, default_symbol="-"):
    if x is None:
        return default_symbol
    if isinstance(x, bool):
        return "True" if x else "False"
    return str(x)

def _print_config_summary(dataset, cfg) -> None:
    """Print all opts, including persisted trim/tximport/align (show defaults if unset)."""
    # Dataset counts
    total = "-"
    try:
        total = len(dataset)
    except Exception:
        pass
    bps = getattr(dataset, "bioprojects", [])
    done = None
    try:
        done = len(dataset.done())
    except Exception:
        done = None

    print("\n============ Run Summary ============")
    print(f"Mode:                 {_val(getattr(dataset, 'mode', None))}")
    print(f"Samples (total):      {_val(total)}")
    print(f"BioProjects (total):  {_val(len(bps))}")
    if done is not None:
        print(f"Samples (done):       {_val(done)}")
    print(f"Output directory:     {_val(getattr(cfg, 'outdir', None))}")
    print(f"Shared directory:     {_val(getattr(cfg, 'shared', None))}")
    print(f"Log file:             {_val(getattr(cfg, 'log', None))}")
    print(f"Threads:              {_val(getattr(cfg, 'threads', None))}")
    print(f"Reference:            {_val(getattr(cfg, 'reference_path', None))}")
    print(f"Tx2Gene:              {_val(getattr(cfg, 'tx2gene', None))}")

    # ---- Trim (fastp) ----
    print(f"Fastp window size:    {_val(getattr(cfg, 'trim_window_size', None))}")
    print(f"Fastp mean quality:   {_val(getattr(cfg, 'trim_mean_quality', None))}")

    # ---- Align (kallisto) ----
    print(f"Align method:         {_val(getattr(cfg, 'align_method', 'kallisto'))}")
    if getattr(cfg, 'align_method', 'kallisto') == "kallisto":
        print(f"Kallisto bootstraps:  {_val(getattr(cfg, 'kallisto_bootstrap', 100))}")
    else:
        print(f"Kallisto bootstraps:  -")

    # ---- tximport ----
    txi_mode = getattr(cfg, "tximport_mode", None)
    txi_ignore = getattr(cfg, "tximport_ignore_tx_version", None)
    print(f"tximport mode:        {_val(txi_mode)}")
    print(f"tximport ignore ver.: {_val(txi_ignore)}")

    print(f"Force re-run:         {_val(getattr(cfg, 'force', False))}")
    print(f"No-exec (dry run):    {_val(getattr(cfg, 'dry_run', False))}")

    use_cache = not getattr(cfg, 'no_cache', False)
    print(f"Use cache:            {_val(use_cache)}")

    # --- Cache volume summary ---
    if use_cache:
        cg = getattr(cfg, "cache_gb", None)
        if cg is None:
            # automatic mode (final high watermark computed by CacheGate)
            # We cannot compute it now (dry run), so show "auto"
            print(f"Cache volume (GiB):   auto")
        else:
            print(f"Cache volume (GiB):   {cg}")
    else:
        print(f"Cache volume (GiB):   -")

    print("=====================================\n")

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
    no_cache: bool,
    cache_gb: int | None = None,
    keep_fastq: bool,
) -> None:
    if input_path is None or reference_path is None:
        raise click.usageError("Missing required options: -i/--input and -r/--reference. See 'hulk -h'.")

    # Load persisted sections from ~/.hulk.json (or fixed path)
    persisted = _cfg_load(None)
    trim_cfg = persisted.get("trim", {}) or {}
    align_cfg = persisted.get("align", {}) or {}
    txi_cfg   = persisted.get("tximport", {}) or {}

    # ---- Trim defaults ----
    DEFAULT_WS = 4
    DEFAULT_MQ = 20
    ws = trim_cfg.get("window_size", DEFAULT_WS)
    mq = trim_cfg.get("mean_quality", DEFAULT_MQ)

    # ---- Align defaults ----
    align_method = (align_cfg.get("method") or "kallisto").lower()
    kallisto_bootstrap = int(align_cfg.get("bootstrap", 100))

    # ---- tximport defaults ----
    txi_mode = txi_cfg.get("mode") or "raw_counts"
    txi_ignore = bool(txi_cfg.get("ignore_tx_version", False))

    # Build Config
    cfg = Config(
        input_path=input_path,
        reference_path=reference_path,
        outdir=output_dir,
        min_threads=min_threads,
        max_threads=max_threads,
        verbose=verbosity,
        force=force,
        aggregate=aggregate,
        dry_run=dry_run,
        tx2gene=tx2gene_path,
        keep_fastq=keep_fastq,

        no_cache=no_cache,
        cache_gb=cache_gb,

        trim_window_size=ws,
        trim_mean_quality=mq,
        align_method=align_method,
        kallisto_bootstrap=kallisto_bootstrap,
        tximport_mode=txi_mode,
        tximport_ignore_tx_version=txi_ignore,
    )

    # Validate reference
    trans_suff = transcriptome_suffixes(reference_path)
    if trans_suff not in TRANSCRIPTOME_FORMATS:
        raise click.UsageError(
            f"Transcriptome file must be one of: {', '.join(sorted(TRANSCRIPTOME_FORMATS))} (got: {reference_path.name})"
        )

    # Dataset
    if input_path.is_dir():
        fastq_files = scan_fastqs(input_path)
        if not fastq_files:
            raise click.UsageError(f"No FASTQ files found in directory: {input_path}")
        dataset = Dataset.from_fastq_dir(input_path, cfg)
        df = None
    else:
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

    # tx2gene validation if provided
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

    # Print ALL opts (incl persisted defaults)
    _print_config_summary(dataset, cfg)

    if dry_run:
        click.secho("\n✅ Dry run complete. No tools executed.\n", fg="green")
        sys.exit(0)

    if not yes:
        click.echo("")
        click.confirm("Proceed with the run? (Y/n)", default=True, abort=True)

    pipeline(dataset, cfg)
    sys.exit(0)

# ---------------------------- Subcommands ----------------------------

@click.group(cls=HulkCommand, help="Run the HULK pipeline with explicit options.")
@click.option("-i", "--input", "input_path",
              type=click.Path(exists=True, dir_okay=True, path_type=Path),
              required=False,
              help="Input table (.csv/.tsv/.txt) with 'Run','BioProject','Model' OR a directory of FASTQ files.")
@click.option("-r", "--reference", "reference_path",
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              required=False,
              help="Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx).")
@click.option("-o", "--output", "output_dir",
              type=click.Path(dir_okay=True, file_okay=False, path_type=Path),
              default=DEFAULT_OUTDIR, show_default=True,
              help="Output directory.")
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
              help="Create a merged TPM table across BioProjects; with --gene-counts, also a global gene-counts table.")
@click.option("-n", "--dry-run", is_flag=True,
              help="Validate inputs and configuration, print plan, and exit without running tools.")
@click.option("-g", "--gene-counts", "tx2gene_path",
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None,
              help="Enable gene counts using a tx2gene (.csv) with columns 'transcript_id','gene_id'.")
@click.option("--keep-fastq", is_flag=True, help="Keep FASTQ files.")
def run_cmd(**kwargs):
    _run_pipeline(**kwargs)

@cli.command("trim", cls=HulkCommand, help="Configure fastp trimming defaults.")
@click.option("-ws", "--window-size", type=int, default=None, help="fastp sliding window size.")
@click.option("-mq", "--mean-quality", type=int, default=None, help="fastp mean quality threshold.")
@click.option("--reset-defaults", is_flag=True, help="Reset trim options to built-in defaults.")
def trim(window_size: int | None, mean_quality: int | None, reset_defaults: bool):
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
@click.option("-m", "--mode",
              type=click.Choice(["raw_counts", "length_scaled_tpm", "scaled_tpm", "dtu_scaled_tpm"], case_sensitive=False),
              default=None,
              help="tximport aggregation/normalization mode.")
@click.option("--ignore-tx-version", "ignore_tx_version", is_flag=True, default=False,
              help="Strip transcript version suffixes before matching.")
@click.option("--reset-defaults", is_flag=True, help="Reset tximport options to built-in defaults.")
def tximport(mode: str | None, ignore_tx_version: bool, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("tximport")
        click.secho(f"Reset tximport settings to defaults at {p}", fg="green")
        return
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

@cli.command("align", cls=HulkCommand, help="Configure alignment/quantification method and options.")
@click.option("--method",
              type=click.Choice(["kallisto"], case_sensitive=False),
              default="kallisto", show_default=True,
              help="Alignment/quantification backend.")
@click.option("-b", "--bootstrap", type=int, default=100, show_default=True,
              help="Number of bootstrap samples for kallisto quant.")
@click.option("--reset-defaults", is_flag=True, help="Reset alignment options to built-in defaults.")
def align(method: str, bootstrap: int, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("align")
        click.secho(f"Reset align settings to defaults at {p}", fg="green")
        return
    payload = {"method": method.lower(), "bootstrap": int(bootstrap)}
    p = _cfg_update("align", payload)
    click.secho(f"Saved align settings to {p}", fg="green")

# ---------------------------- Entry point ----------------------------

def main():
    os.environ.setdefault("COLUMNS", str(WIDE_HELP))
    cli(standalone_mode=True)

if __name__ == "__main__":
    main()
