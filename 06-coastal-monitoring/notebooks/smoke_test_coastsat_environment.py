"""Quick CoastSat environment check for the optional GGOSSS module.

Run with:
    conda activate coastsat
    python notebooks/smoke_test_coastsat_environment.py

This script does not download satellite imagery and does not require FES2022 data.
"""

from pathlib import Path
import sys


def main():
    coastsat_src = Path.home() / "PythonPackages" / "CoastSat-master"
    if coastsat_src.exists():
        sys.path.insert(0, str(coastsat_src))

    from coastsat import SDS_download, SDS_preprocess, SDS_shoreline, SDS_tools  # noqa: F401

    print("coastsat_modules_ok", True)
    print("coastsat_source_exists", coastsat_src.exists())
    print("coastsat_source", coastsat_src)

    try:
        import ee  # noqa: F401

        print("earthengine_api_ok", True)
    except Exception as exc:
        print("earthengine_api_ok", False)
        print("earthengine_api_error", repr(exc))

    try:
        import pyfes  # noqa: F401

        print("pyfes_ok", True)
    except Exception as exc:
        print("pyfes_ok", False)
        print("pyfes_error", repr(exc))

    workspace = Path.home() / "CoastSat_workspace"
    print("workspace_exists", workspace.exists())
    print("workspace", workspace)
    print("roi_folder_exists", (workspace / "ROI").exists())
    print("data_folder_exists", (workspace / "data").exists())


if __name__ == "__main__":
    main()
