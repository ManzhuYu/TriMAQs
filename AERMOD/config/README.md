# Configuration

Copy settings.example.json to settings.json, then set data_root to your dataset location. All path values in the example are relative to that root. A relative data_root is interpreted from the configuration file's folder. The environment variable TRIMAQS_DATA_ROOT can override it.

The paths section gives input/output locations and executables. The settings section contains each stage's years, pollutants and worker count. A key beginning with a stage name belongs to that stage. Keep placeholders such as {year}, {pollutant}, and {grid_ID}: the scripts fill them for each case.

The shared keys source_parameter_folder, boundary, grid_links and compilation_summary connect stages. The export.* paths are used only when preparing the final input dataset. They may point to working results outside the publication data root. Do not commit your local settings.json if it contains personal machine paths.

Use --check to inspect resolved paths without running a stage. Configuration does not supply missing data or infer units; consult the relevant input README before changing an input table.
