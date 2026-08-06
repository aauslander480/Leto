"""ArcPy finalization for component polygons and dissolved management units."""
from __future__ import annotations

from itertools import combinations
import os
from pathlib import Path

import arcpy
import pandas as pd

from .management_units import (
    assign_final_mu_ids,
    summarize_chain_run,
    summarize_mu_chain,
)


COMPATIBILITY_CANDIDATES = (
    "MGMT_CLASS", "MGMT_TYPE", "DB_GROUP", "OWN_CODE", "OWN_TYPE",
    "RIPARIAN", "WETLAND",
)
WEIGHTED_BIOLOGY_FIELDS = (
    "STDAGE_MEAN", "BALIVE_MEAN", "QMD_MEAN", "TPA_MEAN",
    "HOMOGENEITY", "FORTYPE_PCT", "STDAGE_HOM", "BALIVE_HOM",
    "QMD_HOM", "TPA_HOM",
)
QA_FIELD_TYPES = {
    "PART_COUNT": "LONG", "TOTAL_ACRES": "DOUBLE", "IS_MULTIPART": "SHORT",
    "MAX_PAIR_DIST_M": "DOUBLE", "MEAN_PAIR_DIST_M": "DOUBLE",
    "MAX_CENTROID_DIST_M": "DOUBLE", "MU_SPAN_M": "DOUBLE",
    "CHAIN_FLAG": "SHORT", "CHAIN_RATIO": "DOUBLE",
    "SOURCE_SEG_COUNT": "LONG", "MIN_SEG_ACRES": "DOUBLE",
    "MAX_SEG_ACRES": "DOUBLE", "MEAN_SEG_ACRES": "DOUBLE",
}


def finalize_management_outputs(
    component_source,
    component_output,
    unit_output,
    *,
    proximity_enabled,
    maximum_distance_meters,
    warning_chain_ratio,
    qa_folder,
    write_component_pair_csv=False,
    acreage_tolerance_acres=1e-6,
):
    """Create provenance components, routing-safe final IDs, QA, and dissolved MUs."""
    if os.path.normcase(component_source) != os.path.normcase(component_output):
        _delete(component_output)
        arcpy.management.CopyFeatures(component_source, component_output)
    _ensure_component_fields(component_output)

    available = {field.name.upper(): field.name for field in arcpy.ListFields(component_output)}
    compatibility = [available[name] for name in COMPATIBILITY_CANDIDATES if name in available]
    fields = ["OID@", available["SEG_ID"], available["MU_ID"], "Acres", *compatibility]
    records = []
    with arcpy.da.SearchCursor(component_output, fields) as cursor:
        for row in cursor:
            record = dict(zip(fields, row))
            record["SOURCE_OID"] = int(record.pop("OID@"))
            record["SEG_ID"] = int(record.pop(available["SEG_ID"]))
            record["MU_ID"] = int(record.pop(available["MU_ID"]))
            record["BIO_MU_ID"] = record["MU_ID"]
            records.append(record)
    frame = pd.DataFrame(records)
    if frame.empty:
        raise RuntimeError("ManagementSegments_Final contains no polygons.")
    frame = assign_final_mu_ids(
        frame,
        compatibility_fields=compatibility,
        proximity_enabled=proximity_enabled,
    )
    updates = frame.set_index("SOURCE_OID")[["BIO_MU_ID", "MU_ID", "MU_PART"]].to_dict("index")
    with arcpy.da.UpdateCursor(
        component_output, ["OID@", "SEG_ID", "BIO_MU_ID", "MU_ID", "MU_PART", "ParentSEG_ID"]
    ) as cursor:
        for row in cursor:
            values = updates[int(row[0])]
            row[2] = int(values["BIO_MU_ID"])
            row[3] = int(values["MU_ID"])
            row[4] = int(values["MU_PART"])
            row[5] = int(row[1])
            cursor.updateRow(row)

    qa, pair_rows = _calculate_chain_qa(
        component_output,
        maximum_distance_meters=maximum_distance_meters,
        warning_chain_ratio=warning_chain_ratio,
    )
    qa_path = Path(qa_folder) / "Proximity_MU_Chain_QA.csv"
    flags_path = Path(qa_folder) / "Proximity_MU_Chain_Flags.csv"
    qa.to_csv(qa_path, index=False)
    qa.loc[qa["CHAIN_FLAG"].astype(bool)].to_csv(flags_path, index=False)
    pair_path = None
    if write_component_pair_csv:
        pair_path = Path(qa_folder) / "Proximity_MU_Component_Pairs.csv"
        pd.DataFrame(pair_rows).to_csv(pair_path, index=False)

    _delete(unit_output)
    arcpy.management.Dissolve(
        component_output, unit_output, "MU_ID", multi_part="MULTI_PART"
    )
    aggregate = _aggregate_attributes(component_output, compatibility, qa)
    _write_unit_attributes(unit_output, aggregate, component_output, compatibility)

    component_acres = float(frame["Acres"].sum())
    unit_geometry_acres = 0.0
    with arcpy.da.SearchCursor(unit_output, ["SHAPE@"]) as cursor:
        for (geometry,) in cursor:
            unit_geometry_acres += float(geometry.getArea("GEODESIC", "ACRES"))
    discrepancy = unit_geometry_acres - component_acres
    if abs(discrepancy) > acreage_tolerance_acres:
        raise RuntimeError(
            "Management-unit dissolve changed acreage by "
            f"{discrepancy:.9f} acres (tolerance {acreage_tolerance_acres})."
        )

    summary = summarize_chain_run(qa)
    summary.update({
        "component_polygons": int(len(frame)),
        "final_management_units": int(len(aggregate)),
        "component_acres": component_acres,
        "dissolved_geometry_acres": unit_geometry_acres,
        "acreage_discrepancy": discrepancy,
        "chain_qa_csv": str(qa_path),
        "chain_flags_csv": str(flags_path),
        "component_pair_csv": str(pair_path) if pair_path else None,
    })
    return summary


def _ensure_component_fields(feature_class):
    additions = (
        ("BIO_MU_ID", "LONG"), ("MU_PART", "LONG"),
        ("ParentSEG_ID", "LONG"),
    )
    existing = {field.name.upper() for field in arcpy.ListFields(feature_class)}
    for name, kind in additions:
        if name not in existing:
            arcpy.management.AddField(feature_class, name, kind)


def _calculate_chain_qa(feature_class, *, maximum_distance_meters, warning_chain_ratio):
    description = arcpy.Describe(feature_class)
    spatial_reference = description.spatialReference
    if spatial_reference is None or spatial_reference.type != "Projected":
        raise RuntimeError("Chain QA requires a projected coordinate system.")
    meters_per_unit = float(spatial_reference.metersPerUnit)
    groups = {}
    with arcpy.da.SearchCursor(
        feature_class, ["OID@", "MU_ID", "SEG_ID", "Acres", "SHAPE@"]
    ) as cursor:
        for oid, mu_id, seg_id, acres, geometry in cursor:
            groups.setdefault(int(mu_id), []).append({
                "OID": int(oid), "SEG_ID": int(seg_id),
                "Acres": float(acres), "geometry": geometry,
            })

    qa_rows, pair_rows = [], []
    for mu_id in sorted(groups):
        members = groups[mu_id]
        distances = []
        for left, right in combinations(members, 2):
            boundary = float(left["geometry"].distanceTo(right["geometry"])) * meters_per_unit
            centroid = float(
                left["geometry"].trueCentroid.distanceTo(right["geometry"].trueCentroid)
            ) * meters_per_unit
            distances.append({
                "boundary_distance_m": boundary,
                "centroid_distance_m": centroid,
            })
            pair_rows.append({
                "MU_ID": mu_id,
                "LEFT_OID": left["OID"], "RIGHT_OID": right["OID"],
                "LEFT_SEG_ID": left["SEG_ID"], "RIGHT_SEG_ID": right["SEG_ID"],
                "BOUNDARY_DISTANCE_M": boundary,
                "CENTROID_DISTANCE_M": centroid,
            })
        components = pd.DataFrame([
            {"MU_ID": mu_id, "SEG_ID": item["SEG_ID"], "Acres": item["Acres"]}
            for item in members
        ])
        row = summarize_mu_chain(
            components,
            distances,
            maximum_distance_meters=maximum_distance_meters,
            warning_chain_ratio=warning_chain_ratio,
        )
        row.pop("WARNING_CHAIN_FLAG", None)
        qa_rows.append(row)
    return pd.DataFrame(qa_rows), pair_rows


def _aggregate_attributes(feature_class, compatibility, qa):
    available = {field.name.upper(): field.name for field in arcpy.ListFields(feature_class)}
    weighted = [available[name] for name in WEIGHTED_BIOLOGY_FIELDS if name in available]
    dominant = available.get("FORTYPCD_DOM")
    cursor_fields = ["MU_ID", "BIO_MU_ID", "SEG_ID", "Acres", *compatibility, *weighted]
    if dominant and dominant not in cursor_fields:
        cursor_fields.append(dominant)
    rows = []
    with arcpy.da.SearchCursor(feature_class, cursor_fields) as cursor:
        for row in cursor:
            rows.append(dict(zip(cursor_fields, row)))
    frame = pd.DataFrame(rows)
    qa_lookup = qa.set_index("MU_ID").to_dict("index")
    output = {}
    for mu_id, group in frame.groupby("MU_ID", sort=True):
        acres = pd.to_numeric(group["Acres"], errors="raise")
        record = {
            "MU_ID": int(mu_id),
            "BIO_MU_ID": int(group["BIO_MU_ID"].iloc[0]),
            "Acres": float(acres.sum()),
            "SEG_IDS": ",".join(str(v) for v in sorted(group["SEG_ID"].astype(int).unique())),
        }
        for field in compatibility:
            values = group[field].drop_duplicates()
            if len(values) != 1:
                raise RuntimeError(f"MU_ID {mu_id} mixes incompatible {field} values.")
            record[field] = values.iloc[0]
        for field in weighted:
            values = pd.to_numeric(group[field], errors="coerce")
            valid = values.notna()
            record[field] = (
                float((values[valid] * acres[valid]).sum() / acres[valid].sum())
                if valid.any() else None
            )
        if dominant:
            totals = group.assign(_acres=acres).groupby(dominant, dropna=False)["_acres"].sum()
            record[dominant] = sorted(
                totals.items(), key=lambda item: (-item[1], str(item[0]))
            )[0][0]
        record.update(qa_lookup[int(mu_id)])
        output[int(mu_id)] = record
    return output


def _write_unit_attributes(unit_output, aggregate, component_output, compatibility):
    source_types = {field.name: field for field in arcpy.ListFields(component_output)}
    text_fields = {"SEG_IDS": 4000}
    fields_to_add = {
        "BIO_MU_ID": ("LONG", None), "Acres": ("DOUBLE", None),
        **{name: (kind, None) for name, kind in QA_FIELD_TYPES.items()},
        **{name: ("TEXT", length) for name, length in text_fields.items()},
    }
    for field in compatibility:
        source = source_types[field]
        kind = "TEXT" if source.type == "String" else source.type.upper()
        fields_to_add[field] = (kind, source.length if kind == "TEXT" else None)
    for name in WEIGHTED_BIOLOGY_FIELDS + ("FORTYPCD_DOM",):
        if name in source_types:
            source = source_types[name]
            kind = "TEXT" if source.type == "String" else source.type.upper()
            fields_to_add[name] = (kind, source.length if kind == "TEXT" else None)

    existing = {field.name.upper() for field in arcpy.ListFields(unit_output)}
    for name, (kind, length) in fields_to_add.items():
        if name.upper() not in existing:
            kwargs = {"field_length": length} if length else {}
            arcpy.management.AddField(unit_output, name, kind, **kwargs)
    write_fields = [name for name in fields_to_add if name.upper() != "MU_ID"]
    with arcpy.da.UpdateCursor(unit_output, ["MU_ID", *write_fields]) as cursor:
        for row in cursor:
            values = aggregate[int(row[0])]
            for index, field in enumerate(write_fields, start=1):
                value = values.get(field)
                if field in {"IS_MULTIPART", "CHAIN_FLAG"} and value is not None:
                    value = int(bool(value))
                row[index] = value
            cursor.updateRow(row)


def _delete(path):
    if arcpy.Exists(path):
        arcpy.management.Delete(path)
