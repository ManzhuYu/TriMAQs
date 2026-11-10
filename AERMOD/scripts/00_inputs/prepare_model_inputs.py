"""Prepare the shared emission, parameter and supporting input tables."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from trimaqs_config import path
from input_package import prepare_package


def main():
    archive=Path(path('export.source_archive'))
    prepare_package(path('export.source_model_runs'),archive,path('export.destination_inputs'),
        archive/'03_support/boundary/Counties_in_study_area.shp',
        path('export.compilation_summary'),path('export.report_folder'))


if __name__=='__main__':main()
