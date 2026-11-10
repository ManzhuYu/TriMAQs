# Supplementary background concentrations

These optional Appendix C calculations use the final power-plant and NEI industrial emission inventories, source stack/elevation tables, AERMET files, AERMAP receptors and DEMs. They are separate from UOGD concentrations.

Set TRIMAQS_DATA_ROOT. Run run_AERMOD_power_plants.py --category power_plants --workers 20, or run_AERMOD_NEI_point_sources.py --category other_industry --workers 20. Use --prepare-only first to inspect input issues. Full meteorology, DEMs and the model executable must be obtained separately when omitted from the archive.

Background stack/elevation files are in inputs/parameters/background. Files are power_plant_stack_parameters.csv, power_plant_elevation.csv, industrial_stack_parameters.csv and industrial_elevation.csv. Outputs go to outputs/background_model_runs. run_AERMOD_full_calculation.py supplies shared parsing and storage functions; its main baseline-rerun entry point is not needed here.

Power-plant mass defaults to US short tons and its search-radius thresholds use that unit. NEI radii use metric-tonne totals. Fugitive source exclusions and POINTHOR handling follow the final source-type review. This relocated version has not had a full Appendix C numerical rerun. The published UOGD outputs do not include these background concentrations.
