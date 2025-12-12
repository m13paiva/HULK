#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import json
import shutil
from pathlib import Path
from typing import Any, List, Dict, Optional, Tuple

import click
import pandas as pd
from click.core import ParameterSource

from .core import pipeline
from .entities import Config, Dataset
from .post_processing import run_postprocessing
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
    "\n\n"
    "==================================================================================================================================================\n"
    "HULK is a H(igh-volume b)ulk RNA-seq data preprocessing pipeline for NCBI SRA accessions and user-provided FASTQ files.\n"
    "\n"
    "SRA mode (recommended for public data): given an input table (.csv/.tsv/.txt) with columns 'Run', 'BioProject', and 'Model'.\n"
    "\n"
    "FASTQ mode (local FASTQ folders): -i/--input points to a directory where each sample is a subfolder.\n"
    "\n"
    "Configuration subcommands: 'trim', 'align', 'tximport', 'deseq2', and 'plot' control persistence settings.\n"
    "==================================================================================================================================================\n"
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
    def _fmt_rows_with_spacing(
            self,
            formatter: click.HelpFormatter,
            title: str,
            rows: list[tuple[str, str | None]],
    ) -> None:
        rows = [r for r in rows if r]
        if not rows:
            return

        formatter.write(f"{title}:\n")
        for i, (left, right) in enumerate(rows):
            left = left or ""
            right = right or ""
            visible_left_len = len(left)
            gap = RIGHT_COL - visible_left_len
            if gap < 2:
                gap = 2

            right_lines = right.splitlines() or [""]
            first_line = right_lines[0]
            formatter.write(f"  {left}{' ' * gap}{first_line}\n")

            indent = "  " + " " * RIGHT_COL
            for extra in right_lines[1:]:
                formatter.write(f"{indent}{extra}\n")

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
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    required=False,
    help="Input table (.csv/.tsv) or FASTQ directory.",
)
@click.option(
    "-r",
    "--reference",
    "reference_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
    help="Transcriptome FASTA or kallisto index.",
)
@click.option(
    "-o",
    "--output",
    "output_dir",
    type=click.Path(dir_okay=True, file_okay=False, path_type=Path),
    default=DEFAULT_OUTDIR,
    show_default=True,
    help="Output directory.",
)
@click.option("--min-threads", type=int, default=DEFAULT_MIN_THREADS, show_default=True, help="Threads per SRR.")
@click.option("-t", "--max-threads", type=int, default=DEFAULT_MAX_THREADS, show_default=True, help="Total threads.")
@click.option("--verbosity/--no-verbosity", default=True, show_default=True, help="Live progress.")
@click.option("-y", "--yes", is_flag=True, help="Skip prompts.")
@click.option("-f", "--force", is_flag=True, help="Force re-run (overwrite processed data).")
@click.option("-n", "--dry-run", is_flag=True, help="Validate and plan without running.")
@click.option("-g", "--gene-counts", "tx2gene_path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None, help="tx2gene map for gene-level counts.")
@click.option("--no-cache", is_flag=True, help="Disable SRA cache.")
@click.option("-c", "--cache", "cache_gb", type=int, default=None, help="Max SRA cache size (GiB).")
@click.option("--keep-fastq", is_flag=True, help="Keep trimmed FASTQs.")
@click.option(
    "--seq-tech",
    "seq_tech",
    type=str,
    default=None,
    help="Sequencing technology (Required for FASTQ mode).",
)
@click.option(
    "--target-genes",
    "target_genes_files",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    multiple=True,
    help="File(s) containing target genes (one gene per line). Can be used multiple times.",
)
@click.option(
    "--rem-missing-bps",
    is_flag=True,
    default=False,
    help="DANGER: Remove output folders for BioProjects NOT present in the input table.",
)
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
        dry_run: bool,
        tx2gene_path: Path | None,
        no_cache: bool,
        cache_gb: int | None,
        keep_fastq: bool,
        seq_tech: str | None,
        target_genes_files: Tuple[Path],
        rem_missing_bps: bool,
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
        dry_run=dry_run,
        tx2gene_path=tx2gene_path,
        no_cache=no_cache,
        cache_gb=cache_gb,
        keep_fastq=keep_fastq,
        seq_tech=seq_tech,
        target_genes_files=list(target_genes_files) if target_genes_files else None,
        rem_missing_bps=rem_missing_bps,
    )


# ---------------------------- Summary printing ----------------------------
def _print_trailing_newlines(n_bioprojects: int | None = 0) -> None:
    try:
        n = n_bioprojects or 0
        pad_lines = max(4, n + 2)
        click.echo("\n" * pad_lines, nl=False)
    except Exception:
        pass


def _val(x, default_symbol="-"):
    if x is None:
        return default_symbol
    if isinstance(x, bool):
        return "True" if x else "False"
    return str(x)


def _warn(s): return click.style(f"!!! {s} !!!", fg="red", bold=True)


def _print_config_summary(dataset: Dataset, cfg: Config) -> None:
    """Print complete summary of Config and Dataset options."""

    # --- Danger Zone First ---
    if cfg.force or cfg.rem_missing_bps:
        click.echo("\n" + click.style(">>> DANGER ZONE <<<", fg="red", bold=True, blink=True))
        if cfg.force:
            click.echo(_warn(f"FORCE ENABLED: Will overwrite existing data in: {cfg.outdir}"))
        if cfg.rem_missing_bps:
            click.echo(_warn(f"REM-MISSING-BPS ENABLED: Will DELETE any folder in output not in input table!"))
        click.echo("-" * 40)

    print("\n============ Run Summary ============")

    # --- Dataset Stats ---
    total = "-"
    try:
        total = len(dataset)
    except:
        pass

    bps_len = len(getattr(dataset, "bioprojects", []))

    done = None
    try:
        done = len(dataset.done())
    except:
        pass
    mode = _val(getattr(dataset, 'mode', 'Unknown'))
    print(f"Mode:                 {mode}")
    print(f"Samples (total):      {_val(total)}")
    if mode !="SRA":
        print(f"BioProjects (total):  {_val(bps_len)}")
    if done is not None:
        print(f"Samples (done):       {_val(done)}")
    print("-" * 20)

    # --- Paths ---
    print(f"Input:                {_val(cfg.input_path)}")
    print(f"Output directory:     {_val(cfg.outdir)}")
    print(f"Reference:            {_val(cfg.reference_path)}")
    print(f"Tx2Gene:              {_val(cfg.tx2gene)}")

    # --- Scope ---
    if cfg.target_genes_files:
        tgs = [f.name for f in cfg.target_genes_files]
        print(f"Target Genes:         {len(tgs)} files loaded ({', '.join(tgs)})")
    else:
        print(f"Target Genes:         None (Global analysis)")
    print(f"BioProj Filter:       {_val(cfg.bioproject_filter or 'All')}")

    # --- Resources ---
    print(f"Threads:              min={cfg.min_threads}, max={cfg.max_threads}")
    print(f"Cache:                no_cache={cfg.no_cache}, high={cfg.cache_high_gb}G, low={cfg.cache_low_gb}G")

    # --- Tool Configs ---
    print("-" * 20 + " Tools " + "-" * 20)
    print(f"Seq. technology:      {_val(cfg.seq_tech)}")
    print(f"Keep FASTQ:           {_val(cfg.keep_fastq)}")
    print(f"Trim (fastp):         ws={cfg.trim_window_size}, mq={cfg.trim_mean_quality}")
    print(f"Align ({cfg.align_method}):      boot={cfg.kallisto_bootstrap}")
    print(f"tximport:             mode={cfg.tximport_mode}, ignore_ver={cfg.tximport_ignore_tx_version}")

    # --- DESeq2 & Plots ---
    deseq_status = "ENABLED" if cfg.deseq2_vst_enabled else "DISABLED"
    print(f"DESeq2/VST:           {deseq_status} (Var Threshold: {cfg.deseq2_var_threshold})")
    print(f"Plots:                PCA={cfg.plot_pca}, HM={cfg.plot_heatmap}, VarHM={cfg.plot_var_heatmap}")

    print("=====================================\n")


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
        dry_run: bool,
        tx2gene_path: Path | None,
        no_cache: bool,
        cache_gb: int | None = None,
        keep_fastq: bool,
        seq_tech: str | None,
        target_genes_files: List[Path] | None,
        rem_missing_bps: bool,
) -> None:
    if input_path is None or reference_path is None:
        raise click.UsageError("Missing required options: -i/--input and -r/--reference. See 'hulk -h'.")

    # Load persisted sections from ~/.hulk.json (or fixed path)
    persisted = _cfg_load(None)
    trim_cfg = persisted.get("trim", {}) or {}
    align_cfg = persisted.get("align", {}) or {}
    txi_cfg = persisted.get("tximport", {}) or {}
    plot_cfg = persisted.get("plot", {}) or {}

    # ---- Trim defaults ----
    ws = trim_cfg.get("window_size", 4)
    mq = trim_cfg.get("mean_quality", 20)

    # ---- Align defaults ----
    align_method = (align_cfg.get("method") or "kallisto").lower()
    kallisto_bootstrap = int(align_cfg.get("bootstrap", 100))

    # ---- tximport defaults ----
    txi_mode = txi_cfg.get("mode") or "raw_counts"
    txi_ignore = bool(txi_cfg.get("ignore_tx_version", False))

    # ---- Plot defaults ----
    global_pca = bool(plot_cfg.get("global_pca", False))
    global_heatmap = bool(plot_cfg.get("global_heatmap", False))
    global_var_heatmap = bool(plot_cfg.get("global_var_heatmap", False))
    bp_pca = bool(plot_cfg.get("bp_pca", False))
    bp_heatmap = bool(plot_cfg.get("bp_heatmap", False))

    plot_pca = global_pca or bp_pca
    plot_heatmap = global_heatmap or bp_heatmap
    plot_var_heatmap = global_var_heatmap

    # Build Config
    cfg = Config(
        input_path=input_path,
        reference_path=reference_path,
        outdir=output_dir,
        min_threads=min_threads,
        max_threads=max_threads,
        verbose=verbosity,
        force=force,
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

        plot_pca=plot_pca,
        plot_heatmap=plot_heatmap,
        plot_var_heatmap=plot_var_heatmap,

        seq_tech=seq_tech,
        target_genes_files=target_genes_files,
        rem_missing_bps=rem_missing_bps,
    )

    # Validate reference
    trans_suff = transcriptome_suffixes(reference_path)
    if trans_suff not in TRANSCRIPTOME_FORMATS:
        raise click.UsageError(
            f"Transcriptome file must be one of: {', '.join(sorted(TRANSCRIPTOME_FORMATS))} (got: {reference_path.name})"
        )

    # Dataset Loading
    df = None
    if input_path.is_dir():
        # FASTQ mode
        dataset = Dataset.from_fastq_dir(input_path, cfg)
    else:
        # SRA mode
        suf = input_path.suffix.lower()
        if suf not in INPUT_FORMATS:
            raise click.UsageError(f"Input file must be one of: {', '.join(sorted(INPUT_FORMATS))}")
        try:
            if suf == ".csv":
                df = pd.read_csv(input_path, low_memory=False)
            elif suf == ".tsv":
                df = pd.read_csv(input_path, sep="\t", low_memory=False)
            else:
                df = pd.read_table(input_path)
        except Exception as e:
            raise click.ClickException(f"Could not read input file: {e}")

        if not REQUIRED_INPUT_COLS.issubset(df.columns):
            raise click.UsageError(
                "Input file missing columns: " + ", ".join(sorted(REQUIRED_INPUT_COLS.difference(df.columns))))
        df = df[list(REQUIRED_INPUT_COLS)]
        dataset = Dataset.from_dataframe(df, cfg)

    # tx2gene validation
    if tx2gene_path is not None:
        try:
            tx2gene_df = pd.read_csv(tx2gene_path, sep=None, engine="python")
            if not REQUIRED_TX2GENE_COLS.issubset(tx2gene_df.columns):
                raise click.UsageError("tx2gene missing cols: " + ", ".join(REQUIRED_TX2GENE_COLS))
        except Exception as e:
            raise click.ClickException(str(e))

    # FASTQ mode technology check
    if getattr(dataset, "mode", None) == "FASTQ":
        from .align import list_known_seq_techs, _MODEL_PARAMS
        known_techs = list_known_seq_techs()
        seq_tech = getattr(cfg, "seq_tech", None)
        if not seq_tech:
            raise click.UsageError(f"FASTQ mode requires --seq-tech. Known: {', '.join(known_techs)}")
        if seq_tech.strip().upper() not in _MODEL_PARAMS:
            raise click.UsageError(f"Unknown seq-tech '{seq_tech}'. Known: {', '.join(known_techs)}")

    # ------------------------------------------------------------------
    # DANGEROUS LOGIC: rem-missing-bps
    # ------------------------------------------------------------------
    if rem_missing_bps:
        if getattr(dataset, "mode", None) != "SRA":
            click.secho("[WARNING] --rem-missing-bps ignored: Only applicable in SRA mode (input=Table).", fg="yellow")
        else:
            if output_dir.exists():
                expected_bps = set(getattr(dataset, "bioprojects", []))
                existing_items = [p for p in output_dir.iterdir() if p.is_dir()]
                blacklist = {"shared", "fastq_samples"}

                for folder in existing_items:
                    if folder.name in blacklist:
                        continue

                    if folder.name not in expected_bps:
                        msg = f"[DANGER] Removing extraneous BioProject folder: {folder.name}"
                        click.secho(msg, fg="red", bold=True)
                        if not dry_run:
                            try:
                                shutil.rmtree(folder)
                            except Exception as ex:
                                click.secho(f"Failed to remove {folder}: {ex}", fg="red")
            else:
                click.secho("[INFO] Output directory does not exist, nothing to clean.", fg="blue")

    # Print summary (Cleanly passing both objects to the view function)
    _print_config_summary(dataset, cfg)

    if dry_run:
        click.secho("✅ Dry run complete. No tools executed.\n", fg="green")
        sys.exit(0)

    # Extra warning for the reckless
    if (force or rem_missing_bps) and not yes:
        click.secho("WARNING: You have selected DESTRUCTIVE options (Force and/or Rem-Missing).", fg="red", blink=True,
                    bold=True)
        click.confirm("Are you absolutely sure you want to proceed?", default=False, abort=True)
    elif not yes:
        click.confirm("Proceed with the run?", default=True, abort=True)

    pipeline(dataset, cfg)

    if getattr(dataset, "mode", None) == "FASTQ":
        _print_trailing_newlines(0)
    else:
        _print_trailing_newlines(len(getattr(dataset, "bioprojects", None)))
    sys.exit(0)


# ---------------------------- Subcommands ----------------------------

@click.group(cls=HulkCommand, help="Run the HULK pipeline with explicit options.")
def run_cmd(**kwargs):
    pass


@cli.command("trim", cls=HulkCommand, help="Configure fastp trimming defaults.")
@click.option("-ws", "--window-size", type=int, default=None, help="fastp sliding window size.")
@click.option("-mq", "--mean-quality", type=int, default=None, help="fastp mean quality threshold.")
@click.option("--reset-defaults", is_flag=True, help="Reset trim options.")
def trim(window_size: int | None, mean_quality: int | None, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("trim")
        click.secho(f"Reset trim settings to defaults at {p}", fg="green")
        return
    if window_size is None and mean_quality is None:
        click.echo("Try: hulk trim -h")
        return
    p = _cfg_update("trim", {"window_size": window_size, "mean_quality": mean_quality})
    click.secho(f"Saved trim settings to {p}", fg="green")


@cli.command("tximport", cls=HulkCommand, help="Configure tximport aggregation.")
@click.option("-m", "--mode", type=click.Choice(["raw_counts", "length_scaled_tpm", "scaled_tpm", "dtu_scaled_tpm"],
                                                case_sensitive=False), default=None)
@click.option("--ignore-tx-version", "ignore_tx_version", is_flag=True, default=False)
@click.option("--reset-defaults", is_flag=True)
def tximport(mode: str | None, ignore_tx_version: bool, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("tximport")
        click.secho(f"Reset tximport settings at {p}", fg="green")
        return
    ctx = click.get_current_context()
    flag_provided = (ctx.get_parameter_source("ignore_tx_version") == ParameterSource.COMMANDLINE)
    if mode is None and not flag_provided:
        click.echo("Try: hulk tximport -h")
        return
    payload = {}
    if mode is not None: payload["mode"] = mode
    if flag_provided: payload["ignore_tx_version"] = ignore_tx_version
    p = _cfg_update("tximport", payload)
    click.secho(f"Saved tximport settings to {p}", fg="green")


@cli.command("align", cls=HulkCommand, help="Configure alignment/quantification.")
@click.option("--method", type=click.Choice(["kallisto"], case_sensitive=False), default="kallisto", show_default=True)
@click.option("-b", "--bootstrap", type=int, default=100, show_default=True)
@click.option("--reset-defaults", is_flag=True)
def align(method: str, bootstrap: int, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("align")
        click.secho(f"Reset align settings at {p}", fg="green")
        return
    payload = {"method": method.lower(), "bootstrap": int(bootstrap)}
    p = _cfg_update("align", payload)
    click.secho(f"Saved align settings to {p}", fg="green")


@cli.command("deseq2", cls=HulkCommand, help="Configure DESeq2 and variance filtering options.")
@click.option("--enable/--disable", "enabled", default=None, help="Enable or disable DESeq2/VST steps entirely.")
@click.option("--var-threshold", type=float, default=None,
              help="Low variance threshold (default 0.1). Genes below this var are excluded.")
@click.option("--reset-defaults", is_flag=True, help="Reset DESeq2 options to defaults.")
def deseq2(enabled: bool | None, var_threshold: float | None, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("deseq2")
        click.secho(f"Reset DESeq2 settings to defaults at {p}", fg="green")
        return

    if enabled is None and var_threshold is None:
        click.echo("Try: hulk deseq2 -h (to see options)")
        return

    payload: Dict[str, Any] = {}
    if enabled is not None:
        payload["enabled"] = enabled
    if var_threshold is not None:
        payload["var_threshold"] = var_threshold

    p = _cfg_update("deseq2", payload)
    click.secho(f"Saved DESeq2 settings to {p}", fg="green")


@cli.command("plot", cls=HulkCommand, help="Configure plotting behaviour.")
@click.option("--global-pca", type=bool, default=True)
@click.option("--global-heatmap", type=bool, default=True)
@click.option("--global-var-heatmap", type=bool, default=True)
@click.option("--bp-pca", type=bool, default=False)
@click.option("--bp-heatmap", type=bool, default=False)
@click.option("--reset-defaults", is_flag=True)
def plot(global_pca, global_heatmap, global_var_heatmap, bp_pca, bp_heatmap, reset_defaults):
    if reset_defaults:
        p = _cfg_reset("plot")
        click.secho(f"Reset plot settings at {p}", fg="green")
        return
    if all(x is None for x in [global_pca, global_heatmap, global_var_heatmap, bp_pca, bp_heatmap]):
        click.echo("Try: hulk plot -h")
        return
    payload = {}
    if global_pca is not None: payload["global_pca"] = bool(global_pca)
    if global_heatmap is not None: payload["global_heatmap"] = bool(global_heatmap)
    if global_var_heatmap is not None: payload["global_var_heatmap"] = bool(global_var_heatmap)
    if bp_pca is not None: payload["bp_pca"] = bool(bp_pca)
    if bp_heatmap is not None: payload["bp_heatmap"] = bool(bp_heatmap)
    p = _cfg_update("plot", payload)
    click.secho(f"Saved plot settings to {p}", fg="green")


@cli.command("report", cls=HulkCommand, help="Regenerate plots and matrices from existing output.")
@click.option("-o", "--output", "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
              default=DEFAULT_OUTDIR, show_default=True, help="Output directory to scan.")
@click.option("-g", "--gene-counts", "tx2gene_path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None, help="tx2gene map (required).")
@click.option("--target-genes", "target_genes_files", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              multiple=True, help="Target gene list(s) for targeted heatmaps/matrices.")
@click.option("--no-bp-postprocessing", is_flag=True, help="Skip per-BioProject analysis (only run Global).")
def report(output_dir, tx2gene_path, target_genes_files, no_bp_postprocessing):
    """
    Scans an existing output directory and re-runs the post-processing step
    (DESeq2, VST, Plotting, Seidr matrices) without re-running alignment.

    It is SAFE to run this command concurrently with 'hulk run' to generate
    intermediate reports for samples that have finished.
    """
    click.secho("\n[Report] Scanning output directory for finished samples...", fg="yellow")
    click.secho("[Report] Note: It is safe to run this command while the main pipeline is active.", fg="blue")

    # 1. LOAD PERSISTED SETTINGS
    persisted = _cfg_load()
    plot_cfg = persisted.get("plot", {})

    # Only load the supported flags
    global_pca = bool(plot_cfg.get("global_pca", False))
    global_heatmap = bool(plot_cfg.get("global_heatmap", False))
    global_var_heatmap = bool(plot_cfg.get("global_var_heatmap", False))
    bp_pca = bool(plot_cfg.get("bp_pca", False))
    bp_heatmap = bool(plot_cfg.get("bp_heatmap", False))

    # Combine flags for Config
    plot_pca = global_pca or bp_pca
    plot_heatmap = global_heatmap or bp_heatmap
    plot_var_heatmap = global_var_heatmap

    # 2. Initialize Config (WITHOUT the hallucinated args)
    cfg = Config(
        outdir=output_dir,
        tx2gene=tx2gene_path,
        target_genes_files=list(target_genes_files) if target_genes_files else None,
        plots_only_mode=True,
        # Pass the valid settings
        plot_pca=plot_pca,
        plot_heatmap=plot_heatmap,
        plot_var_heatmap=plot_var_heatmap
    )

    # Fallback for tx2gene lookup
    if not cfg.tx2gene:
        saved_tx = persisted.get("general", {}).get("tx2gene")
        if saved_tx:
            cfg.tx2gene = Path(saved_tx).resolve()

    if not cfg.tx2gene:
        click.secho("Error: tx2gene path is missing. Provide it with -g/--gene-counts.", fg="red")
        return

    try:
        # 3. Reconstruct Dataset
        dataset = Dataset.reconstruct_from_output(cfg)

        click.secho(f"[Report] Found {len(dataset)} samples across {len(dataset.bioprojects)} BioProjects.", fg="green")
        click.secho(f"[Report] Mode: {dataset.mode}", fg="green")

        if cfg.target_genes_files:
            tgs = [f.name for f in cfg.target_genes_files]
            click.secho(f"[Report] Using target lists: {', '.join(tgs)}", fg="cyan")

        if no_bp_postprocessing:
            click.secho("[Report] Skipping per-BioProject analysis (--no-bp-postprocessing).", fg="cyan")

        # 4. Run Post-Processing
        # Pass the skip_bp flag correctly
        run_postprocessing(dataset, cfg, skip_bp=no_bp_postprocessing)

        click.secho("\n[Report] Done. Check output/shared/plots for results.", fg="green")

    except FileNotFoundError as e:
        click.secho(f"[Error] {e}", fg="red")
    except Exception as e:
        click.secho(f"[Error] Unexpected failure: {e}", fg="red")

def main():
    os.environ.setdefault("COLUMNS", str(WIDE_HELP))
    cli(standalone_mode=True)


if __name__ == "__main__":
    main()