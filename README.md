# In-silico PyMtCall and Signac comparison

This directory contains a reproducible test workflow that generates synthetic mitochondrial allele counts, runs the Python PyMtCall caller and the R Signac caller on the same input counts, and compares their variant filters and VAFs. The data are simulated for software comparison; they are not biological measurements.

## Files

- `InSillico.py` is the workflow entry point. It generates inputs for each configured depth cap, runs PyMtCall, calls the R script, calculates comparisons, and writes tables and figures.
- `signac_script.R` creates Signac's sparse count matrix, calls `Signac::IdentifyVariants()` and `Signac::AlleleFreq()`, and writes the result tables consumed by Python.
- `README.md` documents assumptions, thresholds, outputs, and statistics.

## Run

Run `InSillico.py` from this directory with Python 3 and R available. The script uses `Rscript` from `PATH`. It needs the Python packages `anndata`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`, and the PyMtCall package. Set the `PYMTCALL_REPOSITORY` environment variable to the local PyMtCall repository directory before running; the script does not contain a machine-specific fallback path. The R installation needs the `Matrix` and `Signac` packages.

The default run uses 1,000 cells, 5% cell dropout, random seed 2026, and depth caps 150, 100, 50, 20, and 10. These settings are constants near the top of `InSillico.py` and can be changed there.

## Synthetic input generation

For each genomic position, the generator chooses two of A/C/G/T to have reads. For each represented base it samples cell-specific mean depths from a lognormal distribution, samples total reads from a Poisson distribution, and clips the result to the selected depth cap plus 30. It splits reads between forward and reverse strands using a beta-binomial draw. A shared randomly selected subset of cells is set to zero at every position to simulate dropout. Reference bases are simulated independently from A/C/G/T.

The outputs are headerless MGATK-style files in `Files/`:

- `mgatk.A.txt`, `mgatk.C.txt`, `mgatk.G.txt`, `mgatk.T.txt`: position, cell barcode, forward count, reverse count.
- `mgatk.coverage.txt`: position, cell barcode, total coverage summed over the four allele files.
- `chrM_refAllele.txt`: position and simulated reference base.

The AnnData passed to PyMtCall also contains a synthetic sparse matrix of 2,000 ATAC peaks. Each cell has cluster-associated marker peaks, shared background peaks, and randomly selected noise peaks. Scanpy computes a UMAP from these counts (`normalize_total`, `log1p`, PCA, neighbors, and UMAP). The same AnnData and cell order are used to display Python and Signac VAFs on identical coordinates. These artificial ATAC features and clusters are only a visualization scaffold: they are generated independently of the simulated mitochondrial variants and do not represent biological associations.

## Caller configurations

### PyMtCall

`process_variants_fast()` runs with:

- minimum reads per strand: 2;
- per-cell coverage threshold: 20.

The minimum-strand requirement determines whether a cell's alternate-allele count contributes to the Python VAF. The coverage threshold sets that cell's Python VAF to zero if site coverage is below 20. PyMtCall's `variant_names` are the variant-level calls it retained from this process.

### Signac

`IdentifyVariants()` runs with `stabilize_variance = TRUE` and `low_coverage_threshold = 10`. A variant passes the Signac high-confidence filter when all of these hold:

- `n_cells_conf_detected >= 5`;
- `strand_correlation >= 0.65`;
- `vmr > 0.01`.

Signac's `AlleleFreq()` calculates alternate reads divided by total reads at the position using the count matrix. It does not apply PyMtCall's per-cell minimum-strand and coverage mask. Consequently, a cell can have positive Signac VAF and zero Python VAF for the same variant; this is a filtering difference, not necessarily a disagreement about the underlying counts.

## Comparison outputs

### Combined run summary: `PyStats/comparison_summary.csv`

One row per depth cap; filtering and VAF comparison fields are combined so the main results can be compared across runs.

Variant-filter fields:

- `candidate_variant_count`: all candidate substitutions evaluated by Signac plus any Python-called variant not present in that set. This is the comparison universe, not a count of called variants.
- `depth_cap`: simulated count cap used to generate this run's data.
- `signac_high_confidence_count`, `python_variant_call_count`: number of variants passing each method's variant-level call criteria.
- `both_pass_count`, `signac_only_count`, `python_only_count`, `neither_pass_count`: the four outcomes for the two binary variant call sets.
- `raw_agreement`: fraction of candidates with matching pass/fail results.
- `cohens_kappa`: pass/fail agreement adjusted for agreement expected from each method's call rate.
- `jaccard_call_set`: intersection divided by union of the two passing variant sets.
- `mcnemar_exact_pvalue`: exact two-sided binomial test of the discordant candidates (`signac_only` versus `python_only`). A small value indicates asymmetric disagreement. It does not test VAF values.
- `python_min_reads_per_strand`, `python_cell_coverage_threshold`, `signac_min_confident_cells`, `signac_min_strand_correlation`, `signac_min_vmr`: thresholds used in that run.

Cell-level VAF/call fields are calculated only for variants that pass both methods' variant-level filters:

- `shared_passing_variant_count`, `matched_cell_count`: number of shared passing variants and cells compared.
- `both_positive_cell_variant_pairs`, `signac_positive_python_zero_pairs`, `python_positive_signac_zero_pairs`, `both_zero_cell_variant_pairs`: totals over variant-by-cell pairs. A positive VAF is treated as a positive cell-level call. These are pair counts, not counts of unique cells across all variants.
- `python_retention_of_signac_positive_pairs`: fraction of Signac-positive variant/cell pairs that also have positive Python VAF.

VAF distribution and test fields:

- `python_median_across_variant_medians`, `signac_median_across_variant_medians`: median across variants of each caller's per-variant median VAF across matched cells.
- `python_sd_across_variant_medians`, `signac_sd_across_variant_medians`: sample standard deviation across those per-variant medians (not cell-level SD or a confidence interval).
- `median_signac_minus_python`: median paired difference between Signac and Python per-variant medians.
- `median_paired_absolute_difference`: median absolute size of those paired differences.
- `spearman_rho_variant_medians`: Spearman rank correlation between the two sets of variant medians; it measures association, not agreement.
- `wilcoxon_statistic`, `wilcoxon_pvalue`: two-sided Wilcoxon signed-rank comparison of the paired per-variant medians. The variants are the paired observations; cells are not treated as independent replicates. A small p-value indicates systematic differences, not their practical size.

The signed-rank test and correlation are descriptive for this simulated set. Variants can share the same genomic position or simulation structure, so they should not be interpreted as independent biological replicates.

### Per-candidate filter results: `PyStats/filter_calls_<depth>.csv`

One row per candidate substitution, with each caller's pass/fail outcome, a four-category result, and Signac's candidate metrics (`n_cells_conf_detected`, `strand_correlation`, and `vmr`). This file is useful for tracing the summary's call-set agreement counts.

### Per-variant VAF results: `PyStats/variant_comparisons_<depth>.csv`

One row per variant passing both callers' variant-level filters. It includes cells compared; cells with positive VAF in both, Signac only, Python only, or neither; per-caller median and SD VAF across cells; median absolute cell-level VAF difference; and cell-level Spearman correlation. `cells_called_by_both` is the requested distinct-cell count for that individual variant. Summing these counts over variants gives cell/variant pairs, not unique cells.

### Comparison figure: `Plots/comparison_<depth>.pdf`

Each depth cap gets one figure with (1) Signac versus Python per-variant median VAFs, (2) the median of those variant medians with plus/minus one SD, and (3) the variant-level pass/fail confusion matrix with Cohen's kappa and the exact McNemar p-value. It is a comparison summary, not a UMAP plot.

### Per-cell VAF maps: `Plots/variant_umaps_<depth>.pdf`

Each depth cap also gets paired UMAP panels for up to six variants that pass both callers' variant-level filters. The variants with the largest median absolute cell-level VAF differences are shown first. Each row shows PyMtCall and Signac values for the same variant, cells, and Scanpy-computed synthetic ATAC embedding; both panels use the same 0-to-1 color scale. These are per-cell VAF overlays, not maps of variant-level call status. The synthetic peak matrix creates display coordinates only, so these plots must not be interpreted as evidence that a variant associates with an ATAC cluster.

### Signac intermediate files: `Signac/count_cap_<depth>/`

- `signac_candidate_stats.csv`: all substitutions evaluated by `IdentifyVariants()` and Signac's statistics for them.
- `signac_high_confidence.csv`: variants passing Signac's high-confidence filters.
- `signac_allele_freq.csv`: cells-by-variants VAF matrix for the high-confidence calls.

The R script pads a single high-confidence variant with one temporary identified variant for `AlleleFreq()` (which requires at least two variants), then removes the padding column from the exported matrix.

## Interpretation cautions

- The methods do not use identical cell-level filters, so VAF distributions are not expected to match exactly.
- The per-variant VAF comparison deliberately restricts to variants passing both variant-level filters; the separate filter comparison evaluates all candidates.
- Python zeros can mean either no positive alternate call after the strand/coverage filters or no reads meeting those filters. They are counted as not called for the cell-level overlap statistics.
- The MGATK input files in `Files/` and comparison outputs from the current workflow are overwritten by a new run. Summary rows are written together to `PyStats/comparison_summary.csv`; per-depth detailed tables and figures are separate. Files left by earlier script versions (such as target-specific stats, UMAP plots, or old filenames in `Signac/`) are not used by the cleaned workflow and are not automatically deleted.
