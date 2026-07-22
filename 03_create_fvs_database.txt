import os
import shutil
import sqlite3
import pandas as pd

# -------------------------------------------------------
# USER INPUTS
# -------------------------------------------------------

#This folder is the one which houses the csv outputs from the last script 'LETO_CSV_PIPELINE'

out_folder = r"D:\forest_condition_2026\FVS\FVS_Inputs"

#Download a blank FVS .db here and add the path to it: https://www.fs.usda.gov/fvs/software/data.shtml

blank_db = r"D:\forest_condition_2026\FVS\FVS_Data.db"

#Name the final output, maybe be worth importing datetime to keep track of these if we run a lot of them

output_db_name = "FVS_Inventory_ManagementUnits_weighted.db"

#These should ideally stay the same

stand_csv_name = "FVS_StandInit.csv"
tree_csv_name = "FVS_TreeInit.csv"

# The following lines affect which version of FVS will run. For example 'SN' = southern FVS. 
# Each state has a unique fvs code as well 

fvs_variant = "SN"
fvs_state_code = 12       # Florida
default_inv_year = 2022
stand_group_name = "All_Stands"

# FVS sampling design for TREE_COUNT already in trees per acre
#These setting affect how FVS calculates the inputed data. Any changes can cause major inaccuracy.
# After much testing I have found these to be correct for what we are doing here.

basal_area_factor = -1
inv_plot_size = 1
brk_dbh = 0
num_plots = 1
nonstk_plots = 0

output_db = os.path.join(out_folder, output_db_name)
stand_csv = os.path.join(out_folder, stand_csv_name)
tree_csv = os.path.join(out_folder, tree_csv_name)

# -------------------------------------------------------
# COPY BLANK DATABASE
# -------------------------------------------------------

if os.path.exists(output_db):
    os.remove(output_db)

shutil.copy2(blank_db, output_db)

print("Copied blank database to:")
print(output_db)

# -------------------------------------------------------
# READ CSVs
# -------------------------------------------------------

stands = pd.read_csv(stand_csv, dtype=str)
trees = pd.read_csv(tree_csv, dtype=str)

# -------------------------------------------------------
# PREPARE STANDS
# -------------------------------------------------------

stand_out = pd.DataFrame()

stand_out["STAND_ID"] = stands["STAND_ID"]
stand_out["STAND_CN"] = stands["STAND_ID"]
stand_out["VARIANT"] = fvs_variant

stand_out["INV_YEAR"] = pd.to_numeric(
    stands["INV_YEAR"],
    errors="coerce"
).fillna(default_inv_year).astype(int)

stand_out["STATE"] = fvs_state_code
stand_out["GROUPS"] = stand_group_name
stand_out["LOCATION"] = stands["MU_ID"].astype(str)

# Initial stand age from FIA dominant-condition STDAGE
if "AGE" not in stands.columns:
    raise ValueError(
        "FVS_StandInit.csv does not contain an AGE column. "
        "Run the updated LETO_CSV_PIPELINE first."
    )

stand_out["AGE"] = pd.to_numeric(
    stands["AGE"],
    errors="coerce"
).round().astype("Int64")

# ------------------------------------------------------------------
# FVS Sampling Design
#
# BASAL_AREA_FACTOR = -1 tells FVS that TREE_COUNT is already
# expressed as Trees Per Acre (TPA). If left blank, FVS will
# re-expand trees, producing greatly inflated stand volumes.
# ------------------------------------------------------------------

stand_out["BASAL_AREA_FACTOR"] = basal_area_factor
stand_out["INV_PLOT_SIZE"] = inv_plot_size
stand_out["BRK_DBH"] = brk_dbh
stand_out["NUM_PLOTS"] = num_plots
stand_out["NONSTK_PLOTS"] = nonstk_plots

# -------------------------------------------------------
# PREPARE TREES
# -------------------------------------------------------

tree_out = pd.DataFrame()

# Stand identifiers
tree_out["STAND_ID"] = trees["STAND_ID"]
tree_out["STAND_CN"] = trees["STAND_ID"]

# Plot identifiers
tree_out["PLOT_ID"] = 1
tree_out["STANDPLOT_ID"] = trees["STAND_ID"].astype(str) + "_1"
tree_out["STANDPLOT_CN"] = tree_out["STANDPLOT_ID"]

# Create sequential tree IDs within each stand
tree_out["TREE_ID"] = (
    trees.groupby("STAND_ID").cumcount() + 1
)

# Tree attributes
tree_out["TREE_COUNT"] = pd.to_numeric(
    trees["TREE_COUNT"],
    errors="coerce"
)

tree_out["SPECIES"] = trees["SPECIES"].astype(str)

tree_out["DIAMETER"] = pd.to_numeric(
    trees["DIAMETER"],
    errors="coerce"
)

tree_out["HT"] = pd.to_numeric(
    trees["HT"],
    errors="coerce"
)

tree_out["CRRATIO"] = pd.to_numeric(
    trees["CRRATIO"],
    errors="coerce"
)

# Remove unusable records
tree_out = tree_out.dropna(subset=[
    "STAND_ID",
    "STAND_CN",
    "STANDPLOT_CN",
    "STANDPLOT_ID",
    "PLOT_ID",
    "TREE_ID",
    "TREE_COUNT",
    "SPECIES",
    "DIAMETER"
]).copy()

# Keep only columns that exist in FVS_TreeInit
tree_out = tree_out[[
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
    "CRRATIO"
]].copy()


# -------------------------------------------------------
# WRITE TO SQLITE DATABASE
# -------------------------------------------------------

conn = sqlite3.connect(output_db)

# Clear existing inventory tables
conn.execute("DELETE FROM FVS_TreeInit")
conn.execute("DELETE FROM FVS_StandInit")
conn.execute("DELETE FROM FVS_PlotInit")

stand_out.to_sql(
    "FVS_StandInit",
    conn,
    if_exists="append",
    index=False
)

tree_out.to_sql(
    "FVS_TreeInit",
    conn,
    if_exists="append",
    index=False
)

conn.commit()

# -------------------------------------------------------
# CHECK COUNTS
# -------------------------------------------------------

stand_count = conn.execute(
    "SELECT COUNT(*) FROM FVS_StandInit"
).fetchone()[0]

tree_count = conn.execute(
    "SELECT COUNT(*) FROM FVS_TreeInit"
).fetchone()[0]

age_check = pd.read_sql_query(
    """
    SELECT
        COUNT(*) AS TOTAL_STANDS,
        COUNT(AGE) AS STANDS_WITH_AGE,
        MIN(AGE) AS MIN_AGE,
        AVG(AGE) AS MEAN_AGE,
        MAX(AGE) AS MAX_AGE
    FROM FVS_StandInit
    """,
    conn
)

conn.close()

print("Database created successfully.")
print("Stand records:", stand_count)
print("Tree records:", tree_count)

print("\nAge values imported to FVS:")
print(age_check.to_string(index=False))

print(output_db)

