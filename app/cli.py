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

# Import Seidr logic
from .seidr import PRESETS, ALGO_MAP
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
    "Configuration subcommands: 'trim', 'align', 'tximport', 'deseq2', 'plot', and 'seidr' control persistence settings.\n"
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
# --- NEW FLAGS ---
@click.option("--no-bp-postprocessing", is_flag=True, help="Skip per-BioProject post-processing.")
@click.option("--no-global-postprocessing", is_flag=True, help="Skip global (all samples) post-processing.")
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
        no_bp_postprocessing: bool,
        no_global_postprocessing: bool,
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
        no_bp_postprocessing=no_bp_postprocessing,
        no_global_postprocessing=no_global_postprocessing,
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
    if mode != "SRA":
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
    print(f"Adv. Plots:           SampleCor={cfg.plot_sample_cor}, Disp={cfg.plot_dispersion}, TopN={cfg.top_n_vars}")

    # --- Seidr ---
    seidr_status = "ENABLED" if getattr(cfg, "seidr_enabled", False) else "DISABLED"
    print(f"Seidr (Network):      {seidr_status}")

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
        no_bp_postprocessing: bool,
        no_global_postprocessing: bool,
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
    # Default to TRUE for globals if not found in JSON
    global_pca = bool(plot_cfg.get("global_pca", True))
    global_heatmap = bool(plot_cfg.get("global_heatmap", True))
    global_var_heatmap = bool(plot_cfg.get("global_var_heatmap", True))
    bp_pca = bool(plot_cfg.get("bp_pca", True))
    bp_heatmap = bool(plot_cfg.get("bp_heatmap", True))
    sample_cor = bool(plot_cfg.get("sample_cor", True))
    dispersion = bool(plot_cfg.get("dispersion", True))
    top_n = int(plot_cfg.get("top_n", 500))

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

        # New Params
        plot_sample_cor=sample_cor,
        plot_dispersion=dispersion,
        top_n_vars=top_n,

        seq_tech=seq_tech,
        target_genes_files=target_genes_files,
        rem_missing_bps=rem_missing_bps,

        # Post-processing control
        no_bp_postprocessing=no_bp_postprocessing,
        no_global_postprocessing=no_global_postprocessing,
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


@cli.command("trim", cls=HulkCommand, help="Configure fastp trimming defaults.")
@click.option("-ws", "--window-size", type=int, default=None, show_default="4",
              help="fastp sliding window size.")
@click.option("-mq", "--mean-quality", type=int, default=None, show_default="20",
              help="fastp mean quality threshold.")
@click.option("--reset-defaults", is_flag=True, help="Reset trim options.")
def trim(window_size: int | None, mean_quality: int | None, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("trim")
        click.secho(f"Reset trim settings to defaults at {p}", fg="green")
        return

    # Load current to merge
    current_cfg = _cfg_load()
    trim_cfg = current_cfg.get("trim", {})

    if window_size is not None: trim_cfg["window_size"] = window_size
    if mean_quality is not None: trim_cfg["mean_quality"] = mean_quality

    p = _cfg_update("trim", trim_cfg)
    click.secho(f"Saved trim settings to {p}", fg="green")

    # Summary
    click.echo("\nCurrent Trim Configuration:")
    click.echo(f"  Window Size:  {trim_cfg.get('window_size', 4)}")
    click.echo(f"  Mean Quality: {trim_cfg.get('mean_quality', 20)}")


@cli.command("tximport", cls=HulkCommand, help="Configure tximport aggregation.")
@click.option("-m", "--mode", type=click.Choice(["raw_counts", "length_scaled_tpm", "scaled_tpm", "dtu_scaled_tpm"],
                                                case_sensitive=False),
              default=None, show_default="raw_counts")
@click.option("--ignore-tx-version/--keep-tx-version", "ignore_tx_version", default=None, show_default="False",
              help="Ignore transcript version (e.g. ENST.1 -> ENST).")
@click.option("--reset-defaults", is_flag=True)
def tximport(mode: str | None, ignore_tx_version: bool | None, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("tximport")
        click.secho(f"Reset tximport settings at {p}", fg="green")
        return

    current_cfg = _cfg_load()
    txi_cfg = current_cfg.get("tximport", {})

    if mode is not None: txi_cfg["mode"] = mode
    if ignore_tx_version is not None: txi_cfg["ignore_tx_version"] = ignore_tx_version

    p = _cfg_update("tximport", txi_cfg)
    click.secho(f"Saved tximport settings to {p}", fg="green")

    # Summary
    click.echo("\nCurrent tximport Configuration:")
    click.echo(f"  Mode:           {txi_cfg.get('mode', 'raw_counts')}")
    click.echo(f"  Ignore Tx Ver:  {txi_cfg.get('ignore_tx_version', False)}")


@cli.command("align", cls=HulkCommand, help="Configure alignment/quantification.")
@click.option("--method", type=click.Choice(["kallisto"], case_sensitive=False),
              default=None, show_default="kallisto")
@click.option("-b", "--bootstrap", type=int, default=None, show_default="100")
@click.option("--reset-defaults", is_flag=True)
def align(method: str | None, bootstrap: int | None, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("align")
        click.secho(f"Reset align settings at {p}", fg="green")
        return

    current_cfg = _cfg_load()
    align_cfg = current_cfg.get("align", {})

    if method is not None: align_cfg["method"] = method.lower()
    if bootstrap is not None: align_cfg["bootstrap"] = int(bootstrap)

    p = _cfg_update("align", align_cfg)
    click.secho(f"Saved align settings to {p}", fg="green")

    # Summary
    click.echo("\nCurrent Alignment Configuration:")
    click.echo(f"  Method:     {align_cfg.get('method', 'kallisto')}")
    click.echo(f"  Bootstraps: {align_cfg.get('bootstrap', 100)}")


@cli.command("deseq2", cls=HulkCommand, help="Configure DESeq2 and variance filtering options.")
@click.option("--enable/--disable", "enabled", default=None, show_default="True",
              help="Enable or disable DESeq2/VST steps entirely.")
@click.option("--var-threshold", type=float, default=None, show_default="0.1",
              help="Low variance threshold. Genes below this var are excluded.")
@click.option("--reset-defaults", is_flag=True, help="Reset DESeq2 options to defaults.")
def deseq2(enabled: bool | None, var_threshold: float | None, reset_defaults: bool):
    if reset_defaults:
        p = _cfg_reset("deseq2")
        click.secho(f"Reset DESeq2 settings to defaults at {p}", fg="green")
        return

    current_cfg = _cfg_load()
    deseq_cfg = current_cfg.get("deseq2", {})

    if enabled is not None: deseq_cfg["enabled"] = enabled
    if var_threshold is not None: deseq_cfg["var_threshold"] = var_threshold

    p = _cfg_update("deseq2", deseq_cfg)
    click.secho(f"Saved DESeq2 settings to {p}", fg="green")

    # Summary
    click.echo("\nCurrent DESeq2 Configuration:")
    click.echo(f"  Enabled:        {deseq_cfg.get('enabled', True)}")
    click.echo(f"  Var Threshold:  {deseq_cfg.get('var_threshold', 0.1)}")


@cli.command("plot", cls=HulkCommand, help="Configure plotting behaviour.")
@click.option("--global-pca", type=bool, default=None, show_default="True")
@click.option("--global-heatmap", type=bool, default=None, show_default="True")
@click.option("--global-var-heatmap", type=bool, default=None, show_default="True")
@click.option("--bp-pca", type=bool, default=None, show_default="False")
@click.option("--bp-heatmap", type=bool, default=None, show_default="False")
@click.option("--sample-cor", type=bool, default=None, show_default="True",
              help="Plot sample-sample correlation heatmap.")
@click.option("--dispersion", type=bool, default=None, show_default="True",
              help="Plot DESeq2 dispersion estimates.")
@click.option("--top-n", type=int, default=None, show_default="500",
              help="Number of top variable genes.")
@click.option("--reset-defaults", is_flag=True)
def plot(global_pca, global_heatmap, global_var_heatmap, bp_pca, bp_heatmap, sample_cor, dispersion, top_n,
         reset_defaults):
    if reset_defaults:
        p = _cfg_reset("plot")
        click.secho(f"Reset plot settings at {p}", fg="green")
        return

    current_cfg = _cfg_load()
    plot_cfg = current_cfg.get("plot", {})

    if global_pca is not None: plot_cfg["global_pca"] = bool(global_pca)
    if global_heatmap is not None: plot_cfg["global_heatmap"] = bool(global_heatmap)
    if global_var_heatmap is not None: plot_cfg["global_var_heatmap"] = bool(global_var_heatmap)
    if bp_pca is not None: plot_cfg["bp_pca"] = bool(bp_pca)
    if bp_heatmap is not None: plot_cfg["bp_heatmap"] = bool(bp_heatmap)
    if sample_cor is not None: plot_cfg["sample_cor"] = bool(sample_cor)
    if dispersion is not None: plot_cfg["dispersion"] = bool(dispersion)
    if top_n is not None: plot_cfg["top_n"] = int(top_n)

    p = _cfg_update("plot", plot_cfg)
    click.secho(f"Saved plot settings to {p}", fg="green")

    # Summary
    click.echo("\nCurrent Plot Configuration:")
    click.echo(f"  Global PCA:     {plot_cfg.get('global_pca', True)}")
    click.echo(f"  Global Heatmap: {plot_cfg.get('global_heatmap', True)}")
    click.echo(f"  Global Var HM:  {plot_cfg.get('global_var_heatmap', True)}")
    click.echo(f"  Sample Cor:     {plot_cfg.get('sample_cor', True)}")
    click.echo(f"  Dispersion:     {plot_cfg.get('dispersion', True)}")
    click.echo(f"  Top N Genes:    {plot_cfg.get('top_n', 500)}")
    click.echo(f"  BP PCA:         {plot_cfg.get('bp_pca', False)}")
    click.echo(f"  BP Heatmap:     {plot_cfg.get('bp_heatmap', False)}")


@cli.command("seidr", cls=HulkCommand, help="Configure Seidr gene network inference settings (persisted).")
@click.option("--enable/--disable", "enabled", default=None, show_default="True",
              help="Enable or disable Seidr analysis in the pipeline.")
@click.option("--preset", type=click.Choice(list(PRESETS.keys()), case_sensitive=False),
              default=None, show_default="BALANCED",
              help="Algorithm preset configuration.")
@click.option("--algo", "algorithms", type=click.Choice(list(ALGO_MAP.keys()), case_sensitive=False),
              multiple=True, help="Manually select specific algorithms (overrides preset).")
@click.option("-b", "--backbone", type=float, default=None, show_default="1.28",
              help="Backbone significance threshold (Fdr).")
@click.option("-w", "--workers", type=int, default=None, show_default="2",
              help="Number of algorithms to run in parallel.")
@click.option("-t", "--target", "targets", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              multiple=True,
              help="Target gene file(s) (persisted default). Pipeline falls back to main --target-genes if this is empty.")
@click.option("--target-mode", type=click.Choice(["both", "main_only", "targeted_only"]),
              default=None, show_default="targeted_only",
              help="Execution scope for targeted analysis.")
@click.option("--reset-defaults", is_flag=True, help="Reset Seidr options to defaults.")
@click.option("--run", is_flag=True, help="Run Seidr analysis immediately using saved settings.")
@click.option("-f", "--force", is_flag=True, help="Force recalculation (ignore cache) if running.")
def seidr(enabled, preset, algorithms, backbone, workers, targets, target_mode, reset_defaults, run, force):
    """
    Configure defaults for the Seidr inference step.
    Use --run to execute the Seidr analysis immediately based on the saved configuration.
    """
    if reset_defaults:
        p = _cfg_reset("seidr")
        click.secho(f"Reset Seidr options to defaults at {p}", fg="green")
        return

    current_cfg = _cfg_load()
    seidr_cfg = current_cfg.get("seidr", {})

    if enabled is not None:
        seidr_cfg["enabled"] = enabled
        click.secho(f"Seidr analysis {'ENABLED' if enabled else 'DISABLED'}.", fg="green" if enabled else "yellow")

    if preset:
        seidr_cfg["preset"] = preset.upper()
        if "algorithms" in seidr_cfg: del seidr_cfg["algorithms"]

    if algorithms:
        seidr_cfg["algorithms"] = [a.upper() for a in algorithms]
        if "preset" in seidr_cfg: del seidr_cfg["preset"]

    if backbone is not None: seidr_cfg["backbone"] = backbone
    if workers is not None: seidr_cfg["workers"] = workers
    if target_mode: seidr_cfg["target_mode"] = target_mode

    if targets:
        seidr_cfg["targets"] = [str(t.resolve()) for t in targets]

    p = _cfg_update("seidr", seidr_cfg)
    click.secho(f"Seidr settings saved to {p}", fg="green")

    # Print summary
    click.echo("\nCurrent Seidr Configuration:")
    click.echo(f"  Enabled:     {seidr_cfg.get('enabled', True)}")
    if "algorithms" in seidr_cfg:
        click.echo(f"  Algorithms:  {', '.join(seidr_cfg['algorithms'])}")
    else:
        click.echo(f"  Preset:      {seidr_cfg.get('preset', 'BALANCED')}")
    click.echo(f"  Backbone:    {seidr_cfg.get('backbone', 1.28)}")
    click.echo(f"  Parallel:    {seidr_cfg.get('workers', 2)} workers")
    saved_targets = seidr_cfg.get("targets", [])
    if saved_targets:
        click.echo(f"  Default Targets: {len(saved_targets)} files")

    # --- EXECUTION BLOCK ---
    if run or force:
        # Load necessary paths from 'general' section to build Config
        general_cfg = current_cfg.get("general", {})
        saved_outdir = general_cfg.get("outdir")
        saved_tx2gene = general_cfg.get("tx2gene")

        if not saved_outdir:
            click.secho("\n[Error] Cannot run Seidr: 'outdir' not found in configuration. Run 'hulk report' first.",
                        fg="red")
            return

        click.secho(f"\n[Seidr] Executing analysis... (Force={force})", fg="magenta", bold=True)

        # Construct simplified Config
        cfg = Config(
            outdir=Path(saved_outdir),
            tx2gene=Path(saved_tx2gene) if saved_tx2gene else None
        )

        try:
            from .seidr import run_seidr
            run_seidr(cfg, force=force)
            click.secho("[Seidr] Done.", fg="green")
        except Exception as e:
            click.secho(f"[Error] Seidr execution failed: {e}", fg="red")

@cli.command("report", cls=HulkCommand,
             help="Regenerate (or generate snapshots while hulk is running) plots and matrices using saved settings.")
@click.option("-o", "--output", "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
              default=DEFAULT_OUTDIR, show_default=True, help="Output directory to scan.")
@click.option("-g", "--gene-counts", "tx2gene_path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              default=None, help="tx2gene map (required).")
@click.option("--target-genes", "target_genes_files", type=click.Path(exists=True, dir_okay=False, path_type=Path),
              multiple=True, help="Target gene list(s) for targeted heatmaps/matrices.")
@click.option("--no-bp-postprocessing", is_flag=True, help="Skip per-BioProject analysis (only run Global).")
@click.option("--no-global-postprocessing", is_flag=True, help="Skip Global analysis (only run BioProjects).")
@click.option("--fast", is_flag=True, help="Skip DESeq2 recalculation (Fast Mode).")
@click.option("-f", "--force", is_flag=True, help="Force full recalculation (overwrites existing matrices).")

def report(output_dir, tx2gene_path, target_genes_files, no_bp_postprocessing, no_global_postprocessing, fast, force):
    """
    Scans output directory and runs post-processing using settings defined in 'hulk plot'.
    """
    click.secho("\n[Report] Scanning output directory...", fg="yellow")

    # 1. LOAD PERSISTED SETTINGS (The Source of Truth)
    persisted = _cfg_load()
    plot_cfg = persisted.get("plot", {})

    # Determine "Fast Mode" logic
    # Force overrides Fast.
    if force:
        plots_only = False
    else:
        plots_only = fast

    # 2. Build Configuration
    cfg = Config(
        outdir=output_dir,
        tx2gene=tx2gene_path,
        target_genes_files=list(target_genes_files) if target_genes_files else None,
        plots_only_mode=plots_only,  # Controlled by Force/Fast flags


        no_bp_postprocessing=no_bp_postprocessing,
        no_global_postprocessing=no_global_postprocessing,

        # Load plot settings directly from persistence
        plot_pca=plot_cfg.get("global_pca", True),
        plot_heatmap=plot_cfg.get("global_heatmap", True),
        plot_var_heatmap=plot_cfg.get("global_var_heatmap", True),
        plot_sample_cor=plot_cfg.get("sample_cor", True),
        plot_dispersion=plot_cfg.get("dispersion", True),
        top_n_vars=plot_cfg.get("top_n", 500)
    )

    # Fallback for tx2gene lookup if not provided in CLI
    if not cfg.tx2gene:
        saved_tx = persisted.get("general", {}).get("tx2gene")
        if saved_tx:
            cfg.tx2gene = Path(saved_tx).resolve()

    if not cfg.tx2gene:
        click.secho("Error: tx2gene path is missing. Provide it with -g/--gene-counts.", fg="red")
        return

    try:
        # 3. Reconstruct Dataset from disk
        dataset = Dataset.reconstruct_from_output(cfg)

        click.secho(f"[Report] Found {len(dataset)} samples across {len(dataset.bioprojects)} BioProjects.", fg="green")
        click.secho(f"[Report] Mode: {dataset.mode}", fg="green")

        if cfg.target_genes_files:
            tgs = [f.name for f in cfg.target_genes_files]
            click.secho(f"[Report] Using target lists: {', '.join(tgs)}", fg="cyan")

        if no_bp_postprocessing:
            click.secho("[Report] Skipping per-BioProject analysis (--no-bp-postprocessing).", fg="cyan")
        if no_global_postprocessing:
            click.secho("[Report] Skipping Global analysis (--no-global-postprocessing).", fg="cyan")

        if force:
            click.secho("[Report] FORCE ENABLED: Recalculating everything (overwriting old files)...", fg="magenta", bold=True)
        elif plots_only:
            click.secho("[Report] Fast mode: Using existing VST matrix (plots only).", fg="yellow")
        else:
            click.secho("[Report] Recalculating expression matrices (this may take a moment)...", fg="magenta")

        # 4. Run Post-Processing
        run_postprocessing(
            dataset,
            cfg,
            skip_bp=no_bp_postprocessing,
            skip_global=no_global_postprocessing
        )

        click.secho("\n[Report] Done. Check output/shared/plots.", fg="green")

    except FileNotFoundError as e:
        click.secho(f"[Error] {e}", fg="red")
    except Exception as e:
        click.secho(f"[Error] Unexpected failure: {e}", fg="red")

def main():
    os.environ.setdefault("COLUMNS", str(WIDE_HELP))
    cli(standalone_mode=True)


if __name__ == "__main__":
    main()