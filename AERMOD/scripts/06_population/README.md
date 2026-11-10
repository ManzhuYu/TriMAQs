# Population datasets

Run build_annual_population.py first, then build_annual_demographics.py. Both need population_spatial.py beside them and the Census_population.zip source archive in inputs/source_data/census. Set TRIMAQS_DATA_ROOT. These scripts allocate source counts, rebuild tract/grid overlaps, apply the study temporal assumptions and write annual nonnegative integer counts. Completed public population files can be used without rebuilding. Twenty workers are configured; lower WORKERS to fit your allocation. Missing income distributions remain unavailable.
