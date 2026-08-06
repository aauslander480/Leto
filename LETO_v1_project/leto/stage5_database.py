import os
import re
import shutil
import sqlite3
from datetime import datetime

import pandas as pd

# -------------------------------------------------------
# USER INPUTS
# -------------------------------------------------------

# Folder containing FVS_StandInit.csv and FVS_TreeInit.csv.
input_folder = r"D:\forest_condition_2026\FVS\FVS_Inputs"

# Parent folder where one uniquely named folder will be created for each run.
runs_root_folder = r"D:\forest_condition_2026\FVS\FVS_Database_Runs"

# Descriptive label included in the run-folder name.
RUN_NAME = "Hard_Ownership_Boundaries"

# Copy the principal CSV and QA inputs into the run folder for reproducibility.
COPY_RUN_INPUTS = True

blank_db = r"D:\forest_condition_2026\FVS\FVS_Data.db"
stand_csv_name = "FVS_StandInit.csv"
tree_csv_name = "FVS_TreeInit.csv"
output_name_prefix = "FVS_Inventory"

fvs_variant = "SN"
fvs_state_code = 12
default_inv_year = 2022
stand_group_name = "All_Stands"

# FVS sampling design. TREE_COUNT is already expressed as trees per acre.
# Populated from the selected regional configuration by pipeline.py.
basal_area_factor = None
inv_plot_size = None
brk_dbh = None
num_plots = None
nonstk_plots = None

# Files copied into Inputs when present.
run_input_names = [
    stand_csv_name,
    tree_csv_name,
    "Segment_PLT_CN_Weights_Filtered.csv",
    "Segment_FVS_Stands_No_Live_Trees.csv",
    "QA_FVS_Run_Summary.txt",
    "QA_Missing_FVS_Species.csv",
    "QA_Missing_FVS_Species_By_Stand.csv",
]

stand_csv = os.path.join(input_folder, stand_csv_name)
tree_csv = os.path.join(input_folder, tree_csv_name)


def safe_filename(value):
    """Return a Windows-safe file or folder name."""
    text = re.sub(r"[^A-Za-z0-9_-]+", "_", str(value).strip())
    return text.strip("_") or "Unknown"


def create_unique_run_folder():
    """Create a timestamped folder that will never overwrite a prior run."""
    os.makedirs(runs_root_folder, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_name = f"{timestamp}_{safe_filename(RUN_NAME)}"
    run_folder = os.path.join(runs_root_folder, base_name)

    counter = 2
    while os.path.exists(run_folder):
        run_folder = os.path.join(
            runs_root_folder,
            f"{base_name}_{counter}",
        )
        counter += 1

    databases_folder = os.path.join(run_folder, "Databases")
    inputs_folder = os.path.join(run_folder, "Inputs")

    os.makedirs(databases_folder)
    os.makedirs(inputs_folder)

    return run_folder, databases_folder, inputs_folder


def require_file(path, label):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} not found:\n{path}")


def prepare_stands(stands):
    required = {"STAND_ID", "INV_YEAR", "MU_ID", "AGE", "DB_GROUP"}
    missing = required - set(stands.columns)
    if missing:
        raise ValueError(
            f"FVS_StandInit.csv is missing fields: {sorted(missing)}"
        )

    out = pd.DataFrame()
    out["STAND_ID"] = stands["STAND_ID"]
    out["STAND_CN"] = stands["STAND_ID"]
    out["VARIANT"] = fvs_variant
    out["INV_YEAR"] = (
        pd.to_numeric(stands["INV_YEAR"], errors="coerce")
        .fillna(default_inv_year)
        .astype(int)
    )
    out["STATE"] = fvs_state_code
    out["GROUPS"] = stand_group_name
    out["LOCATION"] = stands["MU_ID"].astype(str)
    out["AGE"] = (
        pd.to_numeric(stands["AGE"], errors="coerce")
        .round()
        .astype("Int64")
    )

    # Validated sampling design for weighted TREE_COUNT expressed as TPA.
    out["BASAL_AREA_FACTOR"] = basal_area_factor
    out["INV_PLOT_SIZE"] = inv_plot_size
    out["BRK_DBH"] = brk_dbh
    out["NUM_PLOTS"] = num_plots
    out["NONSTK_PLOTS"] = nonstk_plots

    return out


def prepare_trees(trees):
    required_input = {
        "STAND_ID",
        "TREE_COUNT",
        "SPECIES",
        "DIAMETER",
        "HT",
        "CRRATIO",
    }
    missing = required_input - set(trees.columns)
    if missing:
        raise ValueError(
            f"FVS_TreeInit.csv is missing fields: {sorted(missing)}"
        )

    out = pd.DataFrame()
    out["STAND_ID"] = trees["STAND_ID"]
    out["STAND_CN"] = trees["STAND_ID"]
    out["PLOT_ID"] = 1
    out["STANDPLOT_ID"] = trees["STAND_ID"].astype(str) + "_1"
    out["STANDPLOT_CN"] = out["STANDPLOT_ID"]
    out["TREE_ID"] = trees.groupby("STAND_ID").cumcount() + 1
    out["TREE_COUNT"] = pd.to_numeric(
        trees["TREE_COUNT"], errors="coerce"
    )
    out["SPECIES"] = trees["SPECIES"].astype(str)
    out["DIAMETER"] = pd.to_numeric(
        trees["DIAMETER"], errors="coerce"
    )
    out["HT"] = pd.to_numeric(trees["HT"], errors="coerce")
    out["CRRATIO"] = pd.to_numeric(
        trees["CRRATIO"], errors="coerce"
    )

    required_output = [
        "STAND_ID",
        "STAND_CN",
        "STANDPLOT_CN",
        "STANDPLOT_ID",
        "PLOT_ID",
        "TREE_ID",
        "TREE_COUNT",
        "SPECIES",
        "DIAMETER",
    ]
    out = out.dropna(subset=required_output).copy()

    return out[
        [
            "STAND_ID",
            "STAND_CN",
            "STANDPLOT_CN",
            "STANDPLOT_ID",
            "PLOT_ID",
            "TREE_ID",
            "TREE_COUNT",
            "SPECIES",
            "DIAMETER",
            "HT",
            "CRRATIO",
        ]
    ]


def write_database(output_db, stand_out, tree_out):
    shutil.copy2(blank_db, output_db)

    with sqlite3.connect(output_db) as conn:
        conn.execute("DELETE FROM FVS_TreeInit")
        conn.execute("DELETE FROM FVS_StandInit")
        conn.execute("DELETE FROM FVS_PlotInit")

        stand_out.to_sql(
            "FVS_StandInit",
            conn,
            if_exists="append",
            index=False,
        )
        tree_out.to_sql(
            "FVS_TreeInit",
            conn,
            if_exists="append",
            index=False,
        )
        conn.commit()


def copy_run_inputs(inputs_folder):
    copied = []

    if not COPY_RUN_INPUTS:
        return copied

    for name in run_input_names:
        source = os.path.join(input_folder, name)
        if not os.path.exists(source):
            continue

        destination = os.path.join(inputs_folder, name)
        print(f"Copying run input: {name}")
        shutil.copy2(source, destination)
        copied.append(destination)

    return copied


def write_run_settings(
    settings_path,
    run_folder,
    summaries,
    copied_inputs,
    started_at,
):
    finished_at = datetime.now()

    total_stands = sum(item["STANDS"] for item in summaries)
    total_trees = sum(item["TREE_RECORDS"] for item in summaries)
    total_acres = sum(
        item["ACRES"]
        for item in summaries
        if pd.notna(item["ACRES"])
    )

    lines = [
        "LETO FVS DATABASE RUN",
        "=" * 70,
        f"Run name: {RUN_NAME}",
        f"Run folder: {run_folder}",
        f"Started: {started_at:%Y-%m-%d %H:%M:%S}",
        f"Finished: {finished_at:%Y-%m-%d %H:%M:%S}",
        "",
        "FVS SETTINGS",
        "-" * 70,
        f"Variant: {fvs_variant}",
        f"State code: {fvs_state_code}",
        f"Default inventory year: {default_inv_year}",
        f"Stand group: {stand_group_name}",
        "",
        "FVS SAMPLING DESIGN",
        "-" * 70,
        f"BASAL_AREA_FACTOR = {basal_area_factor}",
        f"INV_PLOT_SIZE = {inv_plot_size}",
        f"BRK_DBH = {brk_dbh}",
        f"NUM_PLOTS = {num_plots}",
        f"NONSTK_PLOTS = {nonstk_plots}",
        "",
        "RUN TOTALS",
        "-" * 70,
        f"Database groups: {len(summaries):,}",
        f"Total stands: {total_stands:,}",
        f"Total acres: {total_acres:,.2f}",
        f"Total tree records: {total_trees:,}",
        "",
        "COPIED INPUTS",
        "-" * 70,
    ]

    if copied_inputs:
        lines.extend(copied_inputs)
    else:
        lines.append("None")

    lines.extend(["", "DATABASE GROUPS", "-" * 70])

    for item in summaries:
        lines.append(
            f"{item['DB_GROUP']}: "
            f"{item['STANDS']:,} stands; "
            f"{item['ACRES']:,.2f} acres; "
            f"{item['TREE_RECORDS']:,} tree records"
        )

    with open(settings_path, "w", encoding="utf-8") as file:
        file.write("\n".join(lines) + "\n")


def main():
    started_at = datetime.now()

    require_file(blank_db, "Blank FVS database")
    require_file(stand_csv, "FVS stand CSV")
    require_file(tree_csv, "FVS tree CSV")

    run_folder, databases_folder, inputs_folder = (
        create_unique_run_folder()
    )

    print("=" * 72)
    print("CREATING UNIQUE FVS DATABASE RUN")
    print("=" * 72)
    print(f"Run folder: {run_folder}")

    print("Reading FVS input CSV files...")
    stands = pd.read_csv(stand_csv, dtype=str)
    trees = pd.read_csv(tree_csv, dtype=str)

    stands["DB_GROUP"] = (
        stands["DB_GROUP"]
        .fillna("Unknown")
        .str.strip()
        .replace("", "Unknown")
    )

    summaries = []

    for group_name, group_stands in stands.groupby(
        "DB_GROUP", sort=True
    ):
        ids = set(group_stands["STAND_ID"].astype(str))
        group_trees = trees[
            trees["STAND_ID"].astype(str).isin(ids)
        ].copy()

        stand_out = prepare_stands(group_stands)
        tree_out = prepare_trees(group_trees)

        output_db = os.path.join(
            databases_folder,
            f"{output_name_prefix}_{safe_filename(group_name)}.db",
        )

        write_database(output_db, stand_out, tree_out)

        if "Acres" in group_stands.columns:
            acres = pd.to_numeric(
                group_stands["Acres"], errors="coerce"
            ).sum()
        else:
            acres = float("nan")

        summary = {
            "DB_GROUP": group_name,
            "STANDS": len(stand_out),
            "ACRES": acres,
            "TREE_RECORDS": len(tree_out),
            "OUTPUT_DB": output_db,
        }
        summaries.append(summary)

        print(f"Created: {output_db}")
        print(
            f"  Stands: {len(stand_out):,}; "
            f"acres: {acres:,.2f}; "
            f"trees: {len(tree_out):,}"
        )

    summary_csv = os.path.join(
        run_folder,
        "FVS_Database_Group_Summary.csv",
    )
    pd.DataFrame(summaries).to_csv(summary_csv, index=False)

    copied_inputs = copy_run_inputs(inputs_folder)

    settings_path = os.path.join(run_folder, "Run_Settings.txt")
    write_run_settings(
        settings_path,
        run_folder,
        summaries,
        copied_inputs,
        started_at,
    )

    print("=" * 72)
    print("DATABASE RUN COMPLETE")
    print("=" * 72)
    print(f"Run folder: {run_folder}")
    print(f"Databases: {databases_folder}")
    print(f"Summary: {summary_csv}")
    print(f"Settings: {settings_path}")


if __name__ == "__main__":
    main()