# Validation status

The shared final concentration and source-contribution products were checked separately during HPC processing. The portable repository has syntax checks and focused synthetic source-parameter tests. It has not been executed end to end on the full public archive.

Before running a new full batch, test one case per source category against its archived final concentration output with identical meteorology, receptor geometry, emissions and model executable. Then check annual mass, calendars, receptor uniqueness and category sums. Configuration changes can change numerical results. FASTALL is disabled in the model writers. NO2 uses the archived ARM2/TTRM2 settings; OZONEVAL 40.0 uses the model default µg/m³, not ppb.

The selected other Tier 2 inventory uses the reviewed source list and 2011/2014/2017/2020 anchors. Do not recreate the broader unfiltered inventory. Legacy sensitivity and uncertainty results have not been rerun under the final settings and are not supplied as current validation. CNX comparisons are excluded from this repository.
