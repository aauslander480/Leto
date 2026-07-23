
import os
import sqlite3
import pandas as pd
import numpy as np
import arcpy

# -------------------------------------------------------
# USER INPUTS
# -------------------------------------------------------

# Geodatabase containing ManagementUnits_Final
gdb = r"D:\forest_condition_2026\Stand_delineation\Attempt3\Attempt3.gdb"

# Management-unit feature class
mu_fc = os.path.join(gdb, "ManagementUnits_Final")

# FVS output database downloaded after running FVS
fvs_out_db = r""C:\Users\aca28742\Downloads\bivariate_test.db""

# FIA REF_SPECIES table downloaded from FIA DataMart.
# Required fields:
#   SPCD
#   SPGRPCD
fia_species_csv = (
    r""C:\Users\aca28742\Downloads\FIADB_REFERENCE\REF_SPECIES.csv""
    )
# Output naming
scenario_name = os.path.splitext(
    os.path.basename(fvs_out_db)
)[0]

out_folder = r"D:\forest_condition_2026\FVS\FVS_Outputs"
os.makedirs(out_folder, exist_ok=True)

summary_csv = os.path.join(
    out_folder,
    f"{scenario_name}_Summary_AllYears.csv"
)

classification_csv = os.path.join(
    out_folder,
    f"{scenario_name}_Composition_AllYears.csv"
)

# -------------------------------------------------------
# COMPOSITION SETTINGS
# -------------------------------------------------------

# Composition is based on the proportion of total basal area
# represented by conifer species.

CONIFER_THRESHOLD = 0.70
DECIDUOUS_THRESHOLD = 0.30


# -------------------------------------------------------
# HELPER FUNCTIONS
# -------------------------------------------------------

def normalize_species_code(series):
    """
    Normalize FIA species codes without introducing decimals.
    """
    return (
        pd.to_numeric(series, errors="coerce")
        .astype("Int64")
        .astype("string")
    )


def classify_composition(conifer_fraction):
    """
    Classify a stand from conifer basal-area proportion.
    """
    if pd.isna(conifer_fraction):
        return "Unknown"

    if conifer_fraction >= CONIFER_THRESHOLD:
        return "Coniferous"

    if conifer_fraction <= DECIDUOUS_THRESHOLD:
        return "Deciduous"

    return "Mixed"


def check_required_columns(df, required, table_name):
    """
    Raise an informative error if fields are missing.
    """
    missing = set(required) - set(df.columns)

    if missing:
        raise ValueError(
            f"{table_name} is missing required fields: "
            f"{sorted(missing)}\n"
            f"Available fields: {list(df.columns)}"
        )


# -------------------------------------------------------
# READ FIA SPECIES REFERENCE
# -------------------------------------------------------

species_ref = pd.read_csv(
    fia_species_csv,
    dtype=str,
    low_memory=False
)

check_required_columns(
    species_ref,
    ["SPCD", "SFTWD_HRDWD"],
    "REF_SPECIES.csv"
)

species_ref["SpeciesFIA"] = normalize_species_code(
    species_ref["SPCD"]
)

species_ref["SFTWD_HRDWD"] = (
    species_ref["SFTWD_HRDWD"]
    .astype("string")
    .str.strip()
    .str.upper()
)

species_ref["SPECIES_CLASS"] = species_ref[
    "SFTWD_HRDWD"
].map({
    "S": "Conifer",
    "H": "Deciduous",
    "SOFTWOOD": "Conifer",
    "HARDWOOD": "Deciduous"
})

print("\nSFTWD_HRDWD values in REF_SPECIES:")
print(
    species_ref["SFTWD_HRDWD"]
    .value_counts(dropna=False)
    .head(20)
)

species_lookup = (
    species_ref[
        [
            "SpeciesFIA",
            "SPECIES_CLASS"
        ]
    ]
    .dropna(subset=["SpeciesFIA"])
    .drop_duplicates(subset=["SpeciesFIA"])
)

print(
    "Species codes in FIA lookup:",
    species_lookup["SpeciesFIA"].nunique()
)


# -------------------------------------------------------
# CONNECT TO FVS DATABASE
# -------------------------------------------------------

conn = sqlite3.connect(fvs_out_db)

tables = pd.read_sql_query(
    """
    SELECT name
    FROM sqlite_master
    WHERE type = 'table'
    ORDER BY name
    """,
    conn
)["name"].tolist()

required_tables = {
    "FVS_Summary2",
    "FVS_TreeList"
}

missing_tables = required_tables - set(tables)

if missing_tables:
    conn.close()

    raise ValueError(
        "The FVS output database is missing required tables: "
        f"{sorted(missing_tables)}"
    )


# -------------------------------------------------------
# READ FVS SUMMARY TABLE
# -------------------------------------------------------

summary = pd.read_sql_query(
    """
    SELECT
        StandID,
        Year,
        Age,
        Tpa,
        BA,
        SDI,
        RDSDI,
        CCF,
        TopHt,
        QMD,
        GMD,
        TCuFt,
        MCuFt,
        SCuFt,
        BdFt,
        Acc,
        Mort,
        MAI,
        ForTyp,
        SizeCls,
        StkCls
    FROM FVS_Summary2
    """,
    conn
)

check_required_columns(
    summary,
    [
        "StandID",
        "Year",
        "Age"
    ],
    "FVS_Summary2"
)

summary["StandID"] = summary["StandID"].astype(str)

summary["Year"] = pd.to_numeric(
    summary["Year"],
    errors="coerce"
).astype("Int64")

summary["Age"] = pd.to_numeric(
    summary["Age"],
    errors="coerce"
).astype("Int64")

summary["MU_ID"] = pd.to_numeric(
    summary["StandID"].str.replace(
        "MU_",
        "",
        regex=False
    ),
    errors="coerce"
).astype("Int64")

summary = summary.dropna(
    subset=[
        "StandID",
        "Year",
        "MU_ID"
    ]
).copy()

# Protect against duplicate stand-year records.
duplicate_summary = summary.duplicated(
    subset=[
        "StandID",
        "Year"
    ],
    keep=False
)

if duplicate_summary.any():
    print(
        "WARNING: Duplicate StandID-Year rows in "
        "FVS_Summary2:",
        duplicate_summary.sum()
    )

    summary = summary.drop_duplicates(
        subset=[
            "StandID",
            "Year"
        ],
        keep="first"
    )


# -------------------------------------------------------
# READ FVS SIMULATED TREE LIST
# -------------------------------------------------------

tree_list = pd.read_sql_query(
    """
    SELECT
        StandID,
        Year,
        SpeciesFIA,
        DBH,
        TPA
    FROM FVS_TreeList
    """,
    conn
)

conn.close()

check_required_columns(
    tree_list,
    [
        "StandID",
        "Year",
        "SpeciesFIA",
        "DBH",
        "TPA"
    ],
    "FVS_TreeList"
)

tree_list["StandID"] = tree_list["StandID"].astype(str)

tree_list["Year"] = pd.to_numeric(
    tree_list["Year"],
    errors="coerce"
).astype("Int64")

tree_list["SpeciesFIA"] = normalize_species_code(
    tree_list["SpeciesFIA"]
)

tree_list["DBH"] = pd.to_numeric(
    tree_list["DBH"],
    errors="coerce"
)

tree_list["TPA"] = pd.to_numeric(
    tree_list["TPA"],
    errors="coerce"
)

tree_list = tree_list.dropna(
    subset=[
        "StandID",
        "Year",
        "SpeciesFIA",
        "DBH",
        "TPA"
    ]
).copy()

tree_list = tree_list[
    (tree_list["DBH"] > 0)
    & (tree_list["TPA"] > 0)
].copy()

print("FVS TreeList records:", len(tree_list))

print(
    "TreeList stand-year combinations:",
    tree_list[
        ["StandID", "Year"]
    ].drop_duplicates().shape[0]
)


# -------------------------------------------------------
# CALCULATE TREE-RECORD BASAL AREA
# -------------------------------------------------------

# Basal area per acre:
#
# BA = 0.005454 × DBH² × TPA
#
# DBH is in inches and TPA is trees per acre.

tree_list["TREE_BA"] = (
    0.005454
    * tree_list["DBH"] ** 2
    * tree_list["TPA"]
)


# -------------------------------------------------------
# CLASSIFY SPECIES
# -------------------------------------------------------

tree_list = tree_list.merge(
    species_lookup,
    on="SpeciesFIA",
    how="left",
    validate="many_to_one"
)

unclassified_species = (
    tree_list.loc[
        tree_list["SPECIES_CLASS"].isna(),
        "SpeciesFIA"
    ]
    .dropna()
    .value_counts()
)

print(
    "Tree records without species class:",
    tree_list["SPECIES_CLASS"].isna().sum()
)

if not unclassified_species.empty:
    print("\nMost common unclassified FIA species codes:")
    print(unclassified_species.head(20))


# -------------------------------------------------------
# SUM BASAL AREA BY STAND, YEAR, AND SPECIES CLASS
# -------------------------------------------------------

classified_trees = tree_list[
    tree_list["SPECIES_CLASS"].notna()
].copy()

composition_long = (
    classified_trees
    .groupby(
        [
            "StandID",
            "Year",
            "SPECIES_CLASS"
        ],
        as_index=False
    )
    .agg(
        BASAL_AREA=("TREE_BA", "sum")
    )
)

composition = (
    composition_long
    .pivot_table(
        index=[
            "StandID",
            "Year"
        ],
        columns="SPECIES_CLASS",
        values="BASAL_AREA",
        aggfunc="sum",
        fill_value=0
    )
    .reset_index()
)

composition.columns.name = None

if "Conifer" not in composition.columns:
    composition["Conifer"] = 0.0

if "Deciduous" not in composition.columns:
    composition["Deciduous"] = 0.0

composition = composition.rename(
    columns={
        "Conifer": "CONIFER_BA",
        "Deciduous": "DECIDUOUS_BA"
    }
)

composition["CLASSIFIED_BA"] = (
    composition["CONIFER_BA"]
    + composition["DECIDUOUS_BA"]
)

composition["CONIFER_PCT"] = np.where(
    composition["CLASSIFIED_BA"] > 0,
    composition["CONIFER_BA"]
    / composition["CLASSIFIED_BA"]
    * 100,
    np.nan
)

composition["DECIDUOUS_PCT"] = np.where(
    composition["CLASSIFIED_BA"] > 0,
    composition["DECIDUOUS_BA"]
    / composition["CLASSIFIED_BA"]
    * 100,
    np.nan
)

composition["CONIFER_FRACTION"] = (
    composition["CONIFER_PCT"] / 100
)

composition["COMPOSITION"] = (
    composition["CONIFER_FRACTION"]
    .apply(classify_composition)
)

composition["COMPOSITION_CODE"] = (
    composition["COMPOSITION"]
    .map({
        "Coniferous": 1,
        "Mixed": 2,
        "Deciduous": 3,
        "Unknown": 9
    })
    .astype("Int64")
)


# -------------------------------------------------------
# JOIN COMPOSITION TO FVS SUMMARY
# -------------------------------------------------------

summary = summary.merge(
    composition[
        [
            "StandID",
            "Year",
            "CONIFER_BA",
            "DECIDUOUS_BA",
            "CLASSIFIED_BA",
            "CONIFER_PCT",
            "DECIDUOUS_PCT",
            "COMPOSITION",
            "COMPOSITION_CODE"
        ]
    ],
    on=[
        "StandID",
        "Year"
    ],
    how="left",
    validate="one_to_one"
)

summary["COMPOSITION"] = (
    summary["COMPOSITION"]
    .fillna("Unknown")
)

summary["COMPOSITION_CODE"] = (
    summary["COMPOSITION_CODE"]
    .fillna(9)
    .astype("Int64")
)


# -------------------------------------------------------
# ROUND COMPOSITION FIELDS
# -------------------------------------------------------

for field in [
    "CONIFER_BA",
    "DECIDUOUS_BA",
    "CLASSIFIED_BA",
    "CONIFER_PCT",
    "DECIDUOUS_PCT"
]:
    if field in summary.columns:
        summary[field] = summary[field].round(2)


# -------------------------------------------------------
# EXPORT ALL-YEAR TABLES
# -------------------------------------------------------

summary.to_csv(
    summary_csv,
    index=False
)

classification_fields = [
    "StandID",
    "MU_ID",
    "Year",
    "Age",
    "COMPOSITION",
    "COMPOSITION_CODE",
    "CONIFER_BA",
    "DECIDUOUS_BA",
    "CLASSIFIED_BA",
    "CONIFER_PCT",
    "DECIDUOUS_PCT"
]

summary[
    classification_fields
].to_csv(
    classification_csv,
    index=False
)

print("\nExported full FVS summary:")
print(summary_csv)

print("\nExported composition classification:")
print(classification_csv)

print("\nYears:")
print(sorted(summary["Year"].dropna().unique()))

print("\nRecords:")
print(len(summary))

print("\nComposition totals across all years:")
print(
    summary["COMPOSITION"]
    .value_counts(dropna=False)
)


# -------------------------------------------------------
# CREATE GIS OUTPUT FOR EACH YEAR
# -------------------------------------------------------

join_fields = [
    "Age",
    "COMPOSITION",
    "COMPOSITION_CODE",
    "CONIFER_BA",
    "DECIDUOUS_BA",
    "CLASSIFIED_BA",
    "CONIFER_PCT",
    "DECIDUOUS_PCT",
    "Tpa",
    "BA",
    "SDI",
    "RDSDI",
    "CCF",
    "TopHt",
    "QMD",
    "GMD",
    "TCuFt",
    "MCuFt",
    "SCuFt",
    "BdFt",
    "Acc",
    "Mort",
    "MAI",
    "ForTyp",
    "SizeCls",
    "StkCls"
]

for year in sorted(
    summary["Year"].dropna().unique()
):

    year = int(year)

    print(f"\nProcessing year {year}...")

    year_summary = summary[
        summary["Year"] == year
    ].copy()

    year_csv = os.path.join(
        out_folder,
        f"{scenario_name}_Summary_{year}.csv"
    )

    year_summary.to_csv(
        year_csv,
        index=False
    )

    table_name = arcpy.ValidateTableName(
        f"{scenario_name}_Summary_{year}",
        gdb
    )

    table_gdb = os.path.join(
        gdb,
        table_name
    )

    if arcpy.Exists(table_gdb):
        arcpy.management.Delete(table_gdb)

    arcpy.conversion.TableToTable(
        year_csv,
        gdb,
        table_name
    )

    output_name = arcpy.ValidateTableName(
        f"{scenario_name}_{year}",
        gdb
    )

    out_fc = os.path.join(
        gdb,
        output_name
    )

    if arcpy.Exists(out_fc):
        arcpy.management.Delete(out_fc)

    arcpy.management.CopyFeatures(
        mu_fc,
        out_fc
    )

    arcpy.management.JoinField(
        in_data=out_fc,
        in_field="MU_ID",
        join_table=table_gdb,
        join_field="MU_ID",
        fields=join_fields
    )

    print("Created:", out_fc)

print("\nDone.")

