# Concentration and exposure statistics

Use the final public concentrations and population tables. No model runs are started.

Create a working folder `DATA_ROOT/working/statistics_inputs` containing links or copies of `concentrations` and `source_apportionment` pointing to outputs/concentrations and outputs/source_contributions. These distinct paths are used by the comparison checks.

`python scripts/07_statistics/review_UOGD_statistics_streamed.py --stage concentrations --root DATA_ROOT --processed DATA_ROOT/working/statistics_inputs --workers 2`

Repeat with `--stage exposure`. For source contributions use review_UOGD_results.py with `--stage source_contributions` and the same roots. Tables go to working/statistics_inputs/report_review. Copy completed CSV tables to outputs/report_review for figures. Do not publish checkpoint databases.

These statistics describe modeled concentrations, population weighting and benchmark comparisons; they do not estimate disease risk. Undefined exceedance ratios remain missing.
