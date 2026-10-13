# High-depth truth benchmark

`TruthBenchmark.py` is a separate benchmark workflow. It imports the current `InSillico.py` module and reuses its PyMtCall wrapper, Signac configuration, and variant-name conversion. It does not edit or write into the existing workflow's `Files/`, `Signac/`, `PyStats/`, or `Plots/` directories. Its generated data and outputs live under the ignored `TruthBenchmark/` directory.

## Benchmark design

- 1,000 cells and 1,000 genomic positions.
- Exactly one known non-reference alternate base at each position: 1,000 true variant identities in total. Each true variant is represented across all cells, with cell-specific VAFs. Whole-cell dropout is disabled (`DROPOUT_RATE = 0`) so dropout does not remove true-positive cell/variant pairs from the benchmark.
- The truth dataset is generated once at high depth. Per-cell/site coverage is sampled from a lognormal-Poisson model with arithmetic mean approximately 1,000 reads. The raw high-depth observed alternate fraction is the VAF reference for that cell and variant.
- True cell/variant calls are defined by whether the alternate allele has at least one read in the high-depth truth data. Across all 1,000 positions, the two other possible non-reference substitutions at each position (2,000 in total) are known-negative candidate variants.
- True cell-specific VAF probabilities are sampled from `Beta(1, 1)` (uniform from 0 to 1). Forward/reverse read assignment uses `Beta(10, 10)` strand fractions shared between reference and alternate reads at a cell/position, so strand counts retain a common site/cell bias.
- For each lower-depth run, each high-depth read is independently retained with probability `target depth / 1,000`. Cells are never subsampled; the same 1,000 barcodes remain in every run. Thus target depths are expected per-cell/per-position coverages, and realized coverage varies randomly.
- Depth targets are 150, 100, 50, 20, and 10. The random seed is fixed for reproducibility.

The baseline raw VAF matrix and truth-variant list are written to `TruthBenchmark/high_depth_truth_vaf.npz` and `TruthBenchmark/truth_variants.csv`.

## Caller settings

PyMtCall is invoked using `run_python_caller()` from `InSillico.py`: minimum two alternate reads on each strand and per-cell coverage threshold 20. This preserves the current workflow's Python filtering settings at every target depth, including 20 and 10. At lower targets, this fixed coverage threshold can zero many Python per-cell VAFs; this is an outcome of the setting under test, not an adjusted threshold.

Signac is invoked using the existing `signac_script.R` and `low_coverage_threshold = 10`. Its variant-level truth calls use the same high-confidence criteria as the current workflow: at least five confidently detected cells, strand correlation at least 0.65, and VMR above 0.01.

## Outputs

- `TruthBenchmark/benchmark_summary.csv`: two rows per target depth, one per caller.
- `TruthBenchmark/benchmark_per_variant.csv`: per-true-variant results for each method and depth, including recovery status, VAF errors, and cell-level call counts.
- `TruthBenchmark/accuracy_by_depth.pdf`: cell-call precision, recall, specificity, F1, and mean absolute VAF error against the high-depth truth, each plotted across all target depths for both methods.
- `TruthBenchmark/truth_variants.csv` and `high_depth_truth_vaf.npz`: generated truth labels and high-depth VAF reference values.
- `TruthBenchmark/depth_<n>/input/`: the MGATK-format downsampled data used for that run.
- `TruthBenchmark/depth_<n>/signac/`: Signac intermediate result tables for that run.

### Summary statistics

The benchmark does not report variant-level sensitivity. PyMtCall can retain a variant in its variant list while masking its VAF to zero in cells that fail per-cell filters, so a variant-list pass/fail score would not capture the cell-level behavior being evaluated.

Cell-level confusion metrics are calculated across candidate-variant/cell pairs. For a true variant, its high-depth observed alternate-read status defines truth in each cell; all non-truth substitutions are negative in every cell. `cell_tp`, `cell_fp`, `cell_fn`, `cell_tn`, `cell_precision`, `cell_call_recall`, `cell_specificity`, and `cell_f1` describe the corresponding cell-level positive-VAF calls. Precision, recall, specificity, and F1 are each shown across depths in the plot; higher is better for all four.

VAF metrics use all 1,000 true variants and all 1,000 cells. Missing/filtered method values are treated as zero. `vaf_mae_vs_high_depth_truth` is the mean absolute difference between method VAF and the unthinned high-depth raw VAF; it is shown in the right plot panel, and lower is better. `vaf_rmse` and `vaf_bias_method_minus_truth` give additional error summaries; `vaf_pearson_r` measures overall linear association. VAF error includes cells where a method failed to report a call, so it reflects both read-sampling error and the caller's filtering.


