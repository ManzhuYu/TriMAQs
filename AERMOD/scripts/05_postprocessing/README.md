# Concentrations and source contributions

Run these stages after the required model cases have finished:

```bash
python run.py compile_concentrations --config config/settings.json
python run.py create_concentration_netcdf --config config/settings.json
python run.py create_source_contributions --config config/settings.json
```

`compile_concentrations.py` checks positive-emission source cases and sums their concentrations at common receptors before averaging the selected receptors for each study grid. It writes category-specific daily and annual CSVs under `outputs/compiled`. The shared `compilation_summary` path connects this step to the NetCDF export. Each stage protects existing output files; inspect incomplete results before retrying.

`create_concentration_netcdf.py` combines drilling/completion plus production into Tier_1, compressors plus other Tier 2 facilities into Tier_2, and the two tiers into UOGD. It writes one compressed NetCDF per year in `outputs/concentrations`, with daily and annual NO2, PM2.5, SO2 and VOC. NO2/SO2 are converted from ug/m3 to ppb at the conventional 25 C/1 atm reference; PM2.5/VOC remain ug/m3. Annual values are arithmetic means of calendar-day means. Outside-area cells remain NaN. Known pre-activity inventory absences and production SO2 contribute zero; unexpected missing inputs block processing.

`create_source_contributions.py` divides source concentrations by ML total concentrations for NO2, PM2.5 and SO2. ML inputs use ppb for gases and ug/m3 for PM2.5. No further gas conversion is made. Missing or nonpositive ML denominators produce NaN ratios, and ratios above one are retained. Each yearly file includes daily ratios, ML concentrations, mean daily fractions, paired-day annual ratios and valid-day counts. The two annual statistics are different and should be selected explicitly for analysis.

The ML grid must align with the study coordinates. PM2.5 filenames use `prediction_DATE.nc`; NO2/SO2 use `DATE.nc` in their configured pollutant/1km folders. Optional overlapping PM predictions are compared before processing. The exported datasets are checked after writing. No population-weighted exposure calculations are performed by these three scripts.
