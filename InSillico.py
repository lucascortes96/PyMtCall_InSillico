"""Generate synthetic mtDNA counts and compare PyMtCall with Signac filtering.

This script writes MGATK-format test counts, runs both callers at several depth
caps, and produces filtering/VAF comparison tables. See README.md for the method,
outputs, dependencies, and interpretation of each statistic.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import binomtest, spearmanr, wilcoxon


ROOT = Path(__file__).resolve().parent
PYMTCALL_REPOSITORY = Path(
    os.environ.get(
        "PYMTCALL_REPOSITORY",
        "/Users/lucascortes/Library/CloudStorage/OneDrive-Personal/Documents/Work/Pickett/ATAC-Seq/Py-mTCall",
    )
)
if not PYMTCALL_REPOSITORY.exists():
    raise FileNotFoundError(
        "PyMtCall repository not found. Set PYMTCALL_REPOSITORY to its directory."
    )
sys.path.insert(0, str(PYMTCALL_REPOSITORY))

import Py_mTCall_tools.variant_calling_fast as py_mtcall  # noqa: E402


CELL_COUNT = 1000
DROPOUT_RATE = 0.05
RANDOM_SEED = 2026
DEPTH_CAPS = (150, 100, 50, 20, 10)
PYTHON_MIN_STRAND_COUNT = 2
PYTHON_COVERAGE_THRESHOLD = 20
SIGNAC_LOW_COVERAGE_THRESHOLD = 10
SIGNAC_HIGH_CONFIDENCE_CELLS = 5
SIGNAC_MIN_STRAND_CORRELATION = 0.65
SIGNAC_MIN_VMR = 0.01

BASES = ("A", "C", "G", "T")
ALLELES = ("A", "T", "C", "G")


def generate_synthetic_counts(
    depth_cap: int,
    cell_count: int = CELL_COUNT,
    dropout_rate: float = DROPOUT_RATE,
    seed: int = RANDOM_SEED,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, list[str]]:
    """Simulate allele counts, total coverage, and reference alleles.

    Each genomic position has reads for two randomly selected bases. For each
    represented base, per-cell depths are Poisson sampled around a lognormal
    cell-specific mean, capped at ``depth_cap + 30``; reads are split between
    forward and reverse strands. The same cells drop out at every position.
    """
    rng = np.random.default_rng(seed)
    positions = np.arange(1, cell_count + 1)
    dropout_cells = rng.random(cell_count) < dropout_rate
    barcodes: list[str] = []
    seen_barcodes: set[str] = set()
    while len(barcodes) < cell_count:
        barcode = "".join(rng.choice(BASES, size=16)) + "-1"
        if barcode not in seen_barcodes:
            seen_barcodes.add(barcode)
            barcodes.append(barcode)

    present_bases = np.zeros((cell_count, len(ALLELES)), dtype=bool)
    for position_index in range(cell_count):
        present_bases[position_index, rng.choice(len(ALLELES), size=2, replace=False)] = True

    allele_files: dict[str, pd.DataFrame] = {}
    total_coverage = np.zeros((cell_count, cell_count), dtype=np.int32)
    barcode_array = np.asarray(barcodes)

    for base_index, base in enumerate(ALLELES):
        allele_positions = positions[present_bases[:, base_index]]
        cell_depth_means = rng.lognormal(
            mean=np.log(depth_cap), sigma=2, size=cell_count
        )
        total_reads = rng.poisson(
            lam=cell_depth_means[None, :],
            size=(len(allele_positions), cell_count),
        )
        total_reads = np.clip(total_reads, 0, depth_cap + 30)

        forward_fraction = rng.beta(10, 10, size=total_reads.shape)
        forward_reads = rng.binomial(total_reads, forward_fraction)
        reverse_reads = total_reads - forward_reads
        counts = np.stack((forward_reads, reverse_reads), axis=-1)
        counts[:, dropout_cells, :] = 0
        total_coverage[allele_positions - 1] += counts.sum(axis=-1)

        allele_files[base] = pd.DataFrame(
            {
                0: np.repeat(allele_positions, cell_count),
                1: np.tile(barcode_array, len(allele_positions)),
                2: counts[:, :, 0].ravel(),
                3: counts[:, :, 1].ravel(),
            }
        )

    coverage_file = pd.DataFrame(
        {
            0: np.repeat(positions, cell_count),
            1: np.tile(barcode_array, cell_count),
            2: total_coverage.ravel(),
        }
    )
    reference_file = pd.DataFrame(
        {0: positions, 1: rng.choice(BASES, size=cell_count)}
    )
    return allele_files, coverage_file, reference_file, barcodes


def write_mgatk_inputs(
    allele_files: dict[str, pd.DataFrame],
    coverage_file: pd.DataFrame,
    reference_file: pd.DataFrame,
) -> Path:
    """Write simulated counts in the headerless format expected by both callers."""
    input_dir = ROOT / "Files"
    input_dir.mkdir(parents=True, exist_ok=True)

    for base, dataframe in allele_files.items():
        dataframe.to_csv(input_dir / f"mgatk.{base}.txt", index=False, header=False)
    coverage_file.to_csv(input_dir / "mgatk.coverage.txt", index=False, header=False)
    reference_file.to_csv(input_dir / "chrM_refAllele.txt", index=False, header=False)
    return input_dir


def make_cell_adata(barcodes: list[str], seed: int = RANDOM_SEED) -> ad.AnnData:
    """Create synthetic ATAC peaks and compute a UMAP for variant overlays.

    Cluster-associated marker peaks give the synthetic embedding visible
    structure. These ATAC features are only a visualization scaffold; they are
    independent of the simulated mitochondrial variants.
    """
    rng = np.random.default_rng(seed)
    cell_count = len(barcodes)
    cluster_count = 8
    marker_peaks_per_cluster = 20
    marker_peak_count = cluster_count * marker_peaks_per_cluster
    peak_count = 2_000
    shared_peaks = np.arange(marker_peak_count, marker_peak_count + 100)
    noise_peaks_per_cell = 20
    cluster = rng.integers(0, cluster_count, size=cell_count)

    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[int] = []
    for cell_index, cluster_index in enumerate(cluster):
        marker_peaks = np.arange(
            cluster_index * marker_peaks_per_cluster,
            (cluster_index + 1) * marker_peaks_per_cluster,
        )
        noise_peaks = rng.choice(
            np.arange(marker_peak_count + len(shared_peaks), peak_count),
            size=noise_peaks_per_cell,
            replace=False,
        )
        selected_peaks = np.concatenate((marker_peaks, shared_peaks, noise_peaks))
        row_indices.extend([cell_index] * len(selected_peaks))
        column_indices.extend(selected_peaks.tolist())
        values.extend((rng.poisson(1.5, size=len(selected_peaks)) + 1).tolist())

    peak_matrix = sparse.csr_matrix(
        (values, (row_indices, column_indices)),
        shape=(cell_count, peak_count),
    )
    obs = pd.DataFrame(
        {"synthetic_cluster": pd.Categorical(cluster.astype(str))},
        index=pd.Index(barcodes, name="barcode"),
    )
    adata = ad.AnnData(X=peak_matrix, obs=obs)

    # Compute a genuine Scanpy embedding from the synthetic peak counts.
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    return adata


def run_python_caller(input_dir: Path, cell_adata: ad.AnnData) -> ad.AnnData:
    """Run PyMtCall with the comparison's fixed strand and cell-coverage filters."""
    return py_mtcall.process_variants_fast(
        str(input_dir),
        cell_adata,
        min_strand_count=PYTHON_MIN_STRAND_COUNT,
        coverage_threshold=PYTHON_COVERAGE_THRESHOLD,
    )


def run_signac(input_dir: Path, depth_cap: int) -> Path:
    """Run the companion R script and return its directory of comparison inputs."""
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError("Rscript was not found on PATH; install R and Signac first.")

    output_dir = ROOT / "Signac" / f"count_cap_{depth_cap}"
    output_dir.mkdir(parents=True, exist_ok=True)
    r_script = ROOT / "signac_script.R"
    subprocess.run(
        [
            rscript,
            str(r_script),
            str(input_dir),
            str(SIGNAC_LOW_COVERAGE_THRESHOLD),
            str(output_dir),
        ],
        check=True,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    return output_dir


def py_variant_to_signac_name(variant_name: str) -> str:
    """Convert ``REF-POS-ALT`` PyMtCall labels to Signac's ``POSREF>ALT`` form."""
    reference, position, alternate = variant_name.split("-", maxsplit=2)
    return f"{position}{reference}>{alternate}"


def load_signac_results(
    output_dir: Path, cell_names: pd.Index
) -> tuple[pd.DataFrame, set[str], pd.DataFrame]:
    """Load Signac's candidate stats, high-confidence calls, and cell VAF matrix."""
    candidate_stats = pd.read_csv(output_dir / "signac_candidate_stats.csv")
    high_confidence = pd.read_csv(output_dir / "signac_high_confidence.csv")
    high_confidence_names = set(high_confidence["variant"].astype(str))

    frequency_file = output_dir / "signac_allele_freq.csv"
    signac_vaf = pd.read_csv(frequency_file, index_col="cell")
    signac_vaf.index = signac_vaf.index.astype(str)
    signac_vaf = signac_vaf.reindex(cell_names.astype(str)).fillna(0)
    return candidate_stats, high_confidence_names, signac_vaf


def filtering_statistics(
    candidate_stats: pd.DataFrame,
    python_names: set[str],
    signac_names: set[str],
    depth_cap: int,
) -> tuple[dict[str, int | float], pd.DataFrame]:
    """Compare variant-level call sets and return summary and per-candidate details."""
    candidate_names = set(candidate_stats["variant"].astype(str)) | python_names
    candidates = sorted(candidate_names)
    calls = pd.DataFrame({
        "variant": candidates,
        "signac_high_confidence": [name in signac_names for name in candidates],
        "python_passed_variant_filter": [name in python_names for name in candidates],
    })
    calls["outcome"] = np.select(
        [
            calls["signac_high_confidence"] & calls["python_passed_variant_filter"],
            calls["signac_high_confidence"],
            calls["python_passed_variant_filter"],
        ],
        ["both_pass", "signac_only", "python_only"],
        default="neither_pass",
    )

    candidate_metrics = candidate_stats.drop_duplicates("variant").set_index("variant")
    for metric in ("n_cells_conf_detected", "strand_correlation", "vmr"):
        calls[f"signac_{metric}"] = calls["variant"].map(candidate_metrics[metric])

    both = int((calls["outcome"] == "both_pass").sum())
    signac_only = int((calls["outcome"] == "signac_only").sum())
    python_only = int((calls["outcome"] == "python_only").sum())
    neither = int((calls["outcome"] == "neither_pass").sum())
    total = len(calls)
    signac_pass = both + signac_only
    python_pass = both + python_only
    observed_agreement = (both + neither) / total if total else np.nan
    signac_rate = signac_pass / total if total else np.nan
    python_rate = python_pass / total if total else np.nan
    expected_agreement = (
        signac_rate * python_rate + (1 - signac_rate) * (1 - python_rate)
        if total else np.nan
    )
    discordant = signac_only + python_only

    summary: dict[str, int | float] = {
        "depth_cap": depth_cap,
        "candidate_variant_count": total,
        "signac_high_confidence_count": signac_pass,
        "python_variant_call_count": python_pass,
        "both_pass_count": both,
        "signac_only_count": signac_only,
        "python_only_count": python_only,
        "neither_pass_count": neither,
        "raw_agreement": observed_agreement,
        "cohens_kappa": (
            (observed_agreement - expected_agreement) / (1 - expected_agreement)
            if total and expected_agreement != 1 else np.nan
        ),
        "jaccard_call_set": (
            both / (signac_pass + python_pass - both)
            if signac_pass + python_pass > both else np.nan
        ),
        "mcnemar_exact_pvalue": (
            float(binomtest(signac_only, discordant, 0.5).pvalue)
            if discordant else 1.0
        ),
        "python_min_reads_per_strand": PYTHON_MIN_STRAND_COUNT,
        "python_cell_coverage_threshold": PYTHON_COVERAGE_THRESHOLD,
        "signac_min_confident_cells": SIGNAC_HIGH_CONFIDENCE_CELLS,
        "signac_min_strand_correlation": SIGNAC_MIN_STRAND_CORRELATION,
        "signac_min_vmr": SIGNAC_MIN_VMR,
    }
    return summary, calls


def vaf_statistics(
    adata: ad.AnnData,
    python_names: set[str],
    signac_names: set[str],
    signac_vaf: pd.DataFrame,
    depth_cap: int,
) -> tuple[dict[str, int | float], pd.DataFrame]:
    """Compare per-variant VAF summaries and cell call overlap for shared variants."""
    python_variant_names = list(adata.uns["variant_names"])
    python_index = {
        py_variant_to_signac_name(name): index
        for index, name in enumerate(python_variant_names)
    }
    signac_index = {name: index for index, name in enumerate(signac_vaf.columns)}
    shared = sorted(python_names & signac_names & set(signac_index))
    cell_indices = adata.obs_names.get_indexer(signac_vaf.index)

    variant_rows: list[dict[str, int | float | str]] = []
    totals = {
        "both_positive_cell_variant_pairs": 0,
        "signac_positive_python_zero_pairs": 0,
        "python_positive_signac_zero_pairs": 0,
        "both_zero_cell_variant_pairs": 0,
    }

    for name in shared:
        py_values = np.asarray(
            adata.obsm["variant_vaf"][cell_indices, python_index[name]], dtype=float
        )
        r_values = signac_vaf[name].to_numpy(dtype=float)
        valid = np.isfinite(py_values) & np.isfinite(r_values)
        py_values = py_values[valid]
        r_values = r_values[valid]
        if not len(py_values):
            continue

        py_called = py_values > 0
        r_called = r_values > 0
        both_cells = int(np.sum(py_called & r_called))
        signac_only_cells = int(np.sum(~py_called & r_called))
        python_only_cells = int(np.sum(py_called & ~r_called))
        neither_cells = int(np.sum(~py_called & ~r_called))
        totals["both_positive_cell_variant_pairs"] += both_cells
        totals["signac_positive_python_zero_pairs"] += signac_only_cells
        totals["python_positive_signac_zero_pairs"] += python_only_cells
        totals["both_zero_cell_variant_pairs"] += neither_cells

        rho = np.nan
        if np.unique(py_values).size > 1 and np.unique(r_values).size > 1:
            rho = float(spearmanr(py_values, r_values).statistic)

        variant_rows.append({
            "variant": name,
            "cells_compared": len(py_values),
            "cells_called_by_both": both_cells,
            "cells_called_by_signac_only": signac_only_cells,
            "cells_called_by_python_only": python_only_cells,
            "cells_called_by_neither": neither_cells,
            "python_median_vaf": float(np.median(py_values)),
            "signac_median_vaf": float(np.median(r_values)),
            "python_vaf_sd_across_cells": float(np.std(py_values, ddof=1)) if len(py_values) > 1 else np.nan,
            "signac_vaf_sd_across_cells": float(np.std(r_values, ddof=1)) if len(r_values) > 1 else np.nan,
            "median_absolute_cell_vaf_difference": float(np.median(np.abs(r_values - py_values))),
            "cell_vaf_spearman_rho": rho,
        })

    per_variant = pd.DataFrame(variant_rows)
    n_variants = len(per_variant)
    summary: dict[str, int | float] = {
        "depth_cap": depth_cap,
        "shared_passing_variant_count": n_variants,
        "matched_cell_count": len(signac_vaf),
        **totals,
    }

    signac_positive_pairs = (
        totals["both_positive_cell_variant_pairs"]
        + totals["signac_positive_python_zero_pairs"]
    )
    summary["python_retention_of_signac_positive_pairs"] = (
        totals["both_positive_cell_variant_pairs"] / signac_positive_pairs
        if signac_positive_pairs else np.nan
    )

    if n_variants:
        python_medians = per_variant["python_median_vaf"].to_numpy()
        signac_medians = per_variant["signac_median_vaf"].to_numpy()
        differences = signac_medians - python_medians
        summary.update({
            "python_median_across_variant_medians": float(np.median(python_medians)),
            "python_sd_across_variant_medians": float(np.std(python_medians, ddof=1)) if n_variants > 1 else np.nan,
            "signac_median_across_variant_medians": float(np.median(signac_medians)),
            "signac_sd_across_variant_medians": float(np.std(signac_medians, ddof=1)) if n_variants > 1 else np.nan,
            "median_signac_minus_python": float(np.median(differences)),
            "median_paired_absolute_difference": float(np.median(np.abs(differences))),
            "spearman_rho_variant_medians": (
                float(spearmanr(python_medians, signac_medians).statistic)
                if np.unique(python_medians).size > 1
                and np.unique(signac_medians).size > 1 else np.nan
            ),
        })
        if np.any(differences != 0):
            test = wilcoxon(python_medians, signac_medians, alternative="two-sided")
            summary["wilcoxon_statistic"] = float(test.statistic)
            summary["wilcoxon_pvalue"] = float(test.pvalue)
        else:
            summary["wilcoxon_statistic"] = 0.0
            summary["wilcoxon_pvalue"] = 1.0
    else:
        summary.update({
            "python_median_across_variant_medians": np.nan,
            "python_sd_across_variant_medians": np.nan,
            "signac_median_across_variant_medians": np.nan,
            "signac_sd_across_variant_medians": np.nan,
            "median_signac_minus_python": np.nan,
            "median_paired_absolute_difference": np.nan,
            "spearman_rho_variant_medians": np.nan,
            "wilcoxon_statistic": np.nan,
            "wilcoxon_pvalue": np.nan,
        })
    return summary, per_variant


def save_comparison_plot(
    depth_cap: int,
    filter_summary: dict[str, int | float],
    filter_calls: pd.DataFrame,
    vaf_summary: dict[str, int | float],
    per_variant: pd.DataFrame,
) -> None:
    """Save the requested VAF medians/SD and variant-filter overlap in one figure."""
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    if per_variant.empty:
        axes[0].text(0.5, 0.5, "No variants passed both filters", ha="center", va="center")
        axes[0].set_axis_off()
        axes[1].set_axis_off()
    else:
        python_medians = per_variant["python_median_vaf"].to_numpy()
        signac_medians = per_variant["signac_median_vaf"].to_numpy()
        axes[0].scatter(python_medians, signac_medians, s=12, alpha=0.45)
        axes[0].plot([0, 1], [0, 1], "--", color="gray")
        axes[0].set(
            xlim=(0, 1), ylim=(0, 1),
            xlabel="PyMtCall median VAF per variant",
            ylabel="Signac median VAF per variant",
            title="Matched variant medians",
        )
        axes[0].text(
            0.04, 0.96,
            f"n = {len(per_variant)} variants\n"
            f"Wilcoxon p = {vaf_summary['wilcoxon_pvalue']:.3g}",
            transform=axes[0].transAxes, va="top",
        )

        method_medians = [
            vaf_summary["python_median_across_variant_medians"],
            vaf_summary["signac_median_across_variant_medians"],
        ]
        method_sds = [
            vaf_summary["python_sd_across_variant_medians"],
            vaf_summary["signac_sd_across_variant_medians"],
        ]
        axes[1].errorbar(
            [0, 1], method_medians, yerr=method_sds,
            fmt="o", capsize=6, markersize=8,
        )
        axes[1].set_xticks([0, 1], ["PyMtCall", "Signac"])
        axes[1].set(ylim=(0, 1), ylabel="VAF", title="Median ± SD across variants")
        axes[1].grid(axis="y", alpha=0.25)

    outcome_counts = filter_calls["outcome"].value_counts()
    confusion = np.array([
        [outcome_counts.get("both_pass", 0), outcome_counts.get("signac_only", 0)],
        [outcome_counts.get("python_only", 0), outcome_counts.get("neither_pass", 0)],
    ])
    axes[2].imshow(confusion, cmap="Blues")
    for row in range(2):
        for column in range(2):
            axes[2].text(column, row, str(confusion[row, column]), ha="center", va="center")
    axes[2].set_xticks([0, 1], ["Pass", "Fail"])
    axes[2].set_yticks([0, 1], ["Pass", "Fail"])
    axes[2].set(xlabel="Python filter", ylabel="Signac filter")
    axes[2].set_title(
        f"Variant filter agreement\n"
        f"κ = {filter_summary['cohens_kappa']:.3g}; "
        f"McNemar p = {filter_summary['mcnemar_exact_pvalue']:.3g}"
    )

    figure.suptitle(f"PyMtCall vs Signac comparison (depth cap {depth_cap})")
    figure.tight_layout()
    plot_dir = ROOT / "Plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(plot_dir / f"comparison_{depth_cap}.pdf", bbox_inches="tight")
    plt.close(figure)


def save_variant_umaps(
    adata: ad.AnnData,
    signac_vaf: pd.DataFrame,
    per_variant: pd.DataFrame,
    depth_cap: int,
    max_variants: int = 6,
) -> None:
    """Plot paired caller VAFs on the same UMAP for the most discordant variants."""
    if per_variant.empty:
        print("No variants passed both filters; skipping the VAF UMAP plot.")
        return

    selected = per_variant.sort_values(
        "median_absolute_cell_vaf_difference", ascending=False
    ).head(max_variants)
    embedding = adata.obsm["X_umap"]
    python_names = list(adata.uns["variant_names"])
    python_indices = {
        py_variant_to_signac_name(name): index
        for index, name in enumerate(python_names)
    }

    figure, axes = plt.subplots(
        nrows=len(selected), ncols=2,
        figsize=(10, 3.6 * len(selected)),
        squeeze=False,
        sharex=True,
        sharey=True,
        layout="constrained",
    )
    for row_index, (_, variant_row) in enumerate(selected.iterrows()):
        variant = variant_row["variant"]
        python_values = adata.obsm["variant_vaf"][:, python_indices[variant]]
        signac_values = signac_vaf[variant].to_numpy(dtype=float)

        for column_index, (values, method) in enumerate((
            (python_values, "PyMtCall"),
            (signac_values, "Signac"),
        )):
            axis = axes[row_index, column_index]
            points = axis.scatter(
                embedding[:, 0], embedding[:, 1],
                c=values, cmap="viridis", vmin=0, vmax=1, s=9,
            )
            axis.set_title(f"{variant} — {method}")
            axis.set_xlabel("UMAP 1")
            if column_index == 0:
                axis.set_ylabel("UMAP 2")
        figure.colorbar(
            points,
            ax=axes[row_index].tolist(),
            label="VAF",
            fraction=0.04,
            pad=0.03,
        )

    figure.suptitle(
        f"Per-cell variant VAFs on shared synthetic ATAC UMAP (depth cap {depth_cap})"
    )
    plot_dir = ROOT / "Plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(plot_dir / f"variant_umaps_{depth_cap}.pdf", bbox_inches="tight")
    plt.close(figure)


def compare_one_depth(depth_cap: int) -> dict[str, int | float]:
    """Generate one depth-cap dataset, run both callers, and save its comparisons."""
    print(f"\nDepth cap: {depth_cap}")
    allele_files, coverage, reference, barcodes = generate_synthetic_counts(depth_cap)
    input_dir = write_mgatk_inputs(allele_files, coverage, reference)
    adata = run_python_caller(input_dir, make_cell_adata(barcodes))
    signac_dir = run_signac(input_dir, depth_cap)
    candidate_stats, signac_high_confidence, signac_vaf = load_signac_results(
        signac_dir, adata.obs_names
    )

    python_names = {
        py_variant_to_signac_name(name)
        for name in adata.uns.get("variant_names", [])
    }
    filter_summary, filter_calls = filtering_statistics(
        candidate_stats, python_names, signac_high_confidence, depth_cap
    )
    vaf_summary, per_variant = vaf_statistics(
        adata, python_names, signac_high_confidence, signac_vaf, depth_cap
    )

    output_dir = ROOT / "PyStats"
    output_dir.mkdir(parents=True, exist_ok=True)
    filter_calls.to_csv(output_dir / f"filter_calls_{depth_cap}.csv", index=False)
    per_variant.to_csv(output_dir / f"variant_comparisons_{depth_cap}.csv", index=False)
    save_comparison_plot(depth_cap, filter_summary, filter_calls, vaf_summary, per_variant)
    save_variant_umaps(adata, signac_vaf, per_variant, depth_cap)

    combined_summary = {**filter_summary, **vaf_summary}
    print(
        f"Variant calls: both={filter_summary['both_pass_count']}, "
        f"Signac only={filter_summary['signac_only_count']}, "
        f"Python only={filter_summary['python_only_count']}; "
        f"filter kappa={filter_summary['cohens_kappa']:.4f}."
    )
    print(
        f"VAF: {vaf_summary['shared_passing_variant_count']} shared variants, "
        f"median paired difference={vaf_summary['median_signac_minus_python']:.4f}, "
        f"Wilcoxon p={vaf_summary['wilcoxon_pvalue']:.3g}."
    )
    return combined_summary


def main() -> None:
    """Run all configured depth caps and save one combined summary table."""
    summaries = [compare_one_depth(depth_cap) for depth_cap in DEPTH_CAPS]
    summary_path = ROOT / "PyStats" / "comparison_summary.csv"
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    print(f"\nComparison summary written to {summary_path}")


if __name__ == "__main__":
    main()
