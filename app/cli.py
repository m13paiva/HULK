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
    "\n\n"
    "==================================================================================================================================================\n"
    "HULK is a H(igh-volume b)ulk RNA-seq data preprocessing pipeline for NCBI SRA accessions and user-provided FASTQ files.\n"
    "\n"
    "SRA mode (recommended for public data): given an input table (.csv/.tsv/.txt) with columns 'Run', 'BioProject', and 'Model' (this table is easily\n"
    "obtainable via NCBI's SRA search -> Send to: -> File -> RunInfo. HULK fetches SRR data (prefetch) with safe resume/overwrite, converts to FASTQ\n"
    "(fasterq-dump), performs QC/trimming with fastp, quantifies against a transcriptome using kallisto, optionally aggregates counts/TPM via tximport,\n"
    "and can run DESeq2/VST to produce PCA and expression/variance heatmap plots.\n"
    "\n"
    "FASTQ mode (local FASTQ folders): -i/--input points to a directory where each sample is a subfolder containing 1 (single-end) or 2 (paired-end)\n"
    "FASTQ files. Because in this mode metadata are not available, you MUST provide the sequencing technology via --seq-tech. The value must match \n"
    "one of HULK's known platforms and is used to choose sensible fragment-length defaults for kallisto single-end quantification.\n"
    "\n"
    "Reference transcriptome: -r/--reference accepts a transcriptome FASTA (.fa/.fasta[.gz])or a prebuilt kallisto index (.idx). If a FASTA is\n"
    "provided, HULK builds a shared index once and reuses it across all samples in the run.\n"
    "\n"
    "Outputs: in SRA mode, results are written under <OUTPUT>/<BioProject>/<Run>/; in FASTQ mode, under <OUTPUT>/fastq_samples/<sample_id>/. In both\n"
    "cases a shared log is written to <OUTPUT>/shared/log.txt and, when enabled, a shared cache is used for SRA prefetch data management.\n"
    "\n"
    "Configuration subcommands: 'trim' persists fastp trimming defaults; 'align' configures alignment/kallisto options (e.g. bootstraps); 'tximport'\n"
    "configures aggregation/normalization mode; and 'plot' controls global and per-BioProject PCA/heatmap plotting behaviour.\n"
    "\n"
    "Re-runs are safe: completed samples are skipped unless --force is given; partial or interrupted runs are cleaned up or resumed when possible.\n"
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

            # Support multi-line right-hand text with proper indentation
            right_lines = right.splitlines() or [""]
            first_line = right_lines[0]
            formatter.write(f"  {left}{' ' * gap}{first_line}\n")

            # Subsequent lines: indent to the same column as the first description
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
    help=(
        "SRA mode: input table (.csv/.tsv/.txt) with 'Run','BioProject','Model'.\n"
        "FASTQ mode: directory where each sample is a subfolder containing 1 (single-end)\n"
        "or 2 (paired-end) FASTQ files. Layout (single vs paired) is inferred from\n"
        "the files present."
    ),
)


@click.option(
    "-r",
    "--reference",
    "reference_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
    help="Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx).",
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
@click.option("--min-threads", type=int, default=DEFAULT_MIN_THREADS, show_default=True,
              help="Minimum number of threads per SRR.")
@click.option("-t", "--max-threads", type=int, default=DEFAULT_MAX_THREADS, show_default=True,
              help="Maximum total threads.")
@click.option("--verbosity/--no-verbosity", default=True, show_default=True,
              help="Show live progress bars and console messages (default: on).")
@click.option("-y", "--yes", is_flag=True, help="Assume 'yes' to prompts and run without asking.")
@click.option("-f", "--force", "--overwrite", is_flag=True,
              help="Force re-run: overwrite totally/partially processed SRRs.")
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
@click.option("--keep-fastq", is_flag=True, help="Keep trimmed FASTQ files.")
@click.option(
    "--deseq2/--no-deseq2",
    "deseq2_enabled",
    default=True,
    show_default=True,
    help="Enable DESeq2 normalization + VST expression matrices for downstream plots/network exports.",
)
@click.option(
    "--seq-tech",
    "seq_tech",
    type=str,
    default=None,
    help=(
        "Sequencing technology/platform (e.g. 'illumina'). "
        "Required in FASTQ mode; used to choose kallisto parameters."
    ),
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
    deseq2_enabled: bool,
    seq_tech: str | None,
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
        deseq2_enabled=deseq2_enabled,
        seq_tech=seq_tech,
    )

# ---------------------------- Summary printing ----------------------------
def _print_trailing_newlines(n_bioprojects: int | None = 0) -> None:
    """
    Pad the terminal with blank lines so the shell prompt appears
    after all tqdm progress bars.
    """
    try:
        # One line per BioProject + a couple extra
        n = n_bioprojects or 0
        pad_lines = max(4, n + 2)

        # Print the blank lines without adding an extra newline from click itself
        click.echo("\n" * pad_lines, nl=False)
    except Exception:
        # Never let cosmetic padding break the CLI
        pass

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
    print(f"Seq. technology:      {_val(getattr(cfg, 'seq_tech', None))}")

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

    # ---- DESeq2 / plotting ----
    print(f"DESeq2/VST enabled:   {_val(getattr(cfg, 'deseq2_vst_enabled', None))}")
    pca_flag = getattr(cfg, "plot_pca", None)
    hm_flag = getattr(cfg, "plot_heatmap", None)
    varhm_flag = getattr(cfg, "plot_var_heatmap", None)
    print(f"Plots (PCA/HM/VarHM): {_val(pca_flag)}/{_val(hm_flag)}/{_val(varhm_flag)}")

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
    dry_run: bool,
    tx2gene_path: Path | None,
    no_cache: bool,
    cache_gb: int | None = None,
    keep_fastq: bool,
    deseq2_enabled: bool,
    seq_tech: str | None,
) -> None:
    if input_path is None or reference_path is None:
        raise click.UsageError("Missing required options: -i/--input and -r/--reference. See 'hulk -h'.")

    # Load persisted sections from ~/.hulk.json (or fixed path)
    persisted = _cfg_load(None)
    trim_cfg = persisted.get("trim", {}) or {}
    align_cfg = persisted.get("align", {}) or {}
    txi_cfg   = persisted.get("tximport", {}) or {}
    plot_cfg  = persisted.get("plot", {}) or {}

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

    # ---- Plot defaults (global + per-BioProject) ----
    global_pca = bool(plot_cfg.get("global_pca", False))
    global_heatmap = bool(plot_cfg.get("global_heatmap", False))
    global_var_heatmap = bool(plot_cfg.get("global_var_heatmap", False))
    bp_pca = bool(plot_cfg.get("bp_pca", False))
    bp_heatmap = bool(plot_cfg.get("bp_heatmap", False))

    # For now, Config just needs to know whether these plot types are requested at all.
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

        deseq2_vst_enabled=deseq2_enabled,
        plot_pca=plot_pca,
        plot_heatmap=plot_heatmap,
        plot_var_heatmap=plot_var_heatmap,

        seq_tech=seq_tech,
    )

    # Validate reference
    trans_suff = transcriptome_suffixes(reference_path)
    if trans_suff not in TRANSCRIPTOME_FORMATS:
        raise click.UsageError(
            f"Transcriptome file must be one of: {', '.join(sorted(TRANSCRIPTOME_FORMATS))} (got: {reference_path.name})"
        )

    # Dataset
    if input_path.is_dir():
        # FASTQ mode:
        #   - root directory contains one subdirectory per sample
        #   - each sample directory contains 1 (SE) or 2 (PE) FASTQs
        # Dataset.from_fastq_dir will validate that each subdir has FASTQs.
        dataset = Dataset.from_fastq_dir(input_path, cfg)
        df = None
    else:
        suf = input_path.suffix.lower()
        if suf not in INPUT_FORMATS:
            raise click.UsageError(
                f"Input file must be one of: {', '.join(sorted(INPUT_FORMATS))} "
                f"(got: {input_path.name})"
            )
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

    # ------------------------------------------------------------------
    # FASTQ mode: require a known sequencing technology
    # ------------------------------------------------------------------
    if getattr(dataset, "mode", None) == "FASTQ":
        from .align import list_known_seq_techs, _MODEL_PARAMS

        known_techs = list_known_seq_techs()
        seq_tech = getattr(cfg, "seq_tech", None)

        if not seq_tech:
            # No --seq-tech given → print all known techs and abort
            click.secho(
                "[WARNING] No sequencing technology provided (--seq-tech).\n"
                "It is required in FASTQ mode to correctly set fragment length parameters "
                "for kallisto quantification.\n",
                fg="yellow",
            )
            click.echo("Known technologies:\n  - " + "\n  - ".join(known_techs))
            raise click.UsageError(
                "Missing --seq-tech. Please specify one of the known sequencing technologies listed above."
            )

        seq_tech_upper = seq_tech.strip().upper()
        if seq_tech_upper not in _MODEL_PARAMS:
            # Provided value not in known presets → print all known techs and abort
            click.secho(
                f"[WARNING] Sequencing technology '{seq_tech}' not recognised in preset table.\n"
                "HULK only accepts known technologies in FASTQ mode.\n",
                fg="yellow",
            )
            click.echo("Known technologies:\n  - " + "\n  - ".join(known_techs))
            raise click.UsageError(
                f"Unknown sequencing technology '{seq_tech}'. "
                "Please use one of the known technologies listed above."
            )


    # Print ALL opts (incl persisted defaults)
    _print_config_summary(dataset, cfg)

    if not cfg.deseq2_vst_enabled and (cfg.plot_pca or cfg.plot_heatmap or cfg.plot_var_heatmap):
        click.secho("Note: DESeq2 is disabled; any configured expression plots will be skipped downstream.", fg="yellow")

    if dry_run:
        click.secho("\n✅ Dry run complete. No tools executed.\n", fg="green")
        sys.exit(0)

    if not yes:
        click.echo("")
        click.confirm("Proceed with the run? (Y/n)", default=True, abort=True)

    pipeline(dataset, cfg)
    if getattr(dataset, "mode", None) == "FASTQ":
        _print_trailing_newlines(0)
    else:
        _print_trailing_newlines(len(getattr(dataset, "bioprojects",None)))

    sys.exit(0)

# ---------------------------- Subcommands ----------------------------

@click.group(cls=HulkCommand, help="Run the HULK pipeline with explicit options.")
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
    required=False,
    help=(
        "SRA mode: input table (.csv/.tsv/.txt) with 'Run','BioProject','Model'.\n"
        "FASTQ mode: directory where each sample is a subfolder containing 1 (single-end)\n"
        "or 2 (paired-end) FASTQ files. Layout (single vs paired) is inferred from\n"
        "the files present."
    ),
)
@click.option(
    "-r",
    "--reference",
    "reference_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=False,
    help="Reference transcriptome (.fasta/.fa[.gz]) or kallisto index (.idx).",
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
@click.option("--min-threads", type=int, default=DEFAULT_MIN_THREADS, show_default=True,
              help="Minimum number of threads per SRR.")
@click.option("-t", "--max-threads", type=int, default=DEFAULT_MAX_THREADS, show_default=True,
              help="Maximum total threads.")
@click.option("--verbosity/--no-verbosity", default=True, show_default=True,
              help="Show live progress bars and console messages (default: on).")
@click.option("-y", "--yes", is_flag=True, help="Assume 'yes' to prompts and run without asking.")
@click.option("-f", "--force", "--overwrite", is_flag=True,
              help="Force re-run: overwrite totally/partially processed SRRs.")
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
@click.option(
    "--deseq2/--no-deseq2",
    "deseq2_enabled",
    default=True,
    show_default=True,
    help="Enable DESeq2 normalization + VST expression matrices for downstream plots/network exports.",
)
@click.option(
    "--seq-tech",
    "seq_tech",
    type=str,
    default=None,
    help=(
        "Sequencing technology/platform (e.g. 'illumina'). "
        "Required in FASTQ mode; used to choose kallisto parameters."
    ),
)
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

@cli.command("plot", cls=HulkCommand, help="Configure global and per-BioProject plotting behaviour.")
@click.option(
    "--global-pca",
    type=bool,
    default=True,
    metavar="BOOL",
    help="Enable/disable global PCA plot across all samples (true/false).",
)
@click.option(
    "--global-heatmap",
    type=bool,
    default=True,
    metavar="BOOL",
    help="Enable/disable global expression heatmap across all samples (true/false).",
)
@click.option(
    "--global-var-heatmap",
    type=bool,
    default=True,
    metavar="BOOL",
    help="Enable/disable global variance heatmap (gene x BioProject) (true/false).",
)
@click.option(
    "--bp-pca",
    type=bool,
    default=False,
    metavar="BOOL",
    help="Enable/disable per-BioProject PCA plots (one PCA per BioProject) (true/false).",
)
@click.option(
    "--bp-heatmap",
    type=bool,
    default=False,
    metavar="BOOL",
    help="Enable/disable per-BioProject expression heatmaps (one heatmap per BioProject) (true/false).",
)
@click.option("--reset-defaults", is_flag=True, help="Reset plot options to built-in defaults.")
def plot(
    global_pca: bool | None,
    global_heatmap: bool | None,
    global_var_heatmap: bool | None,
    bp_pca: bool | None,
    bp_heatmap: bool | None,
    reset_defaults: bool,
):
    """
    Persistently configure which plots are requested.

    Note: plotting only takes effect when DESeq2/VST is enabled for a run
    (via --deseq2/--no-deseq2 on the main command or 'run' subcommand).
    """
    if reset_defaults:
        p = _cfg_reset("plot")
        click.secho(f"Reset plot settings to defaults at {p}", fg="green")
        return

    if (
        global_pca is None
        and global_heatmap is None
        and global_var_heatmap is None
        and bp_pca is None
        and bp_heatmap is None
    ):
        click.echo("Try: hulk plot -h   (to see options)")
        return

    payload: Dict[str, Any] = {}
    if global_pca is not None:
        payload["global_pca"] = bool(global_pca)
    if global_heatmap is not None:
        payload["global_heatmap"] = bool(global_heatmap)
    if global_var_heatmap is not None:
        payload["global_var_heatmap"] = bool(global_var_heatmap)
    if bp_pca is not None:
        payload["bp_pca"] = bool(bp_pca)
    if bp_heatmap is not None:
        payload["bp_heatmap"] = bool(bp_heatmap)

    p = _cfg_update("plot", payload)
    click.secho(f"Saved plot settings to {p}", fg="green")

# ---------------------------- Entry point ----------------------------

def main():
    os.environ.setdefault("COLUMNS", str(WIDE_HELP))
    cli(standalone_mode=True)

if __name__ == "__main__":
    main()
