\# LETO



LETO is a spatial forest landscape modeling workflow that integrates

TreeMap, Forest Inventory and Analysis data, ArcGIS Pro, and the Forest

Vegetation Simulator.



\## Workflow



1\. Delineate landscape management units.

2\. Assign area-weighted FIA donor plots using TreeMap.

3\. Construct synthetic FVS inventories.

4\. Run Forest Vegetation Simulator scenarios.

5\. Join projected forest conditions back to GIS polygons.



\## Main scripts



| Script | Purpose |

|---|---|

| `01\_delineate\_stands.py` | Creates management units and TreeMap donor weights |

| `02\_create\_fvs\_inputs.py` | Creates weighted FVS stand and tree inventories |

| `03\_create\_fvs\_database.py` | Loads the inventory into an FVS SQLite database |

| `04\_join\_fvs\_outputs.py` | Converts FVS output into GIS-ready layers |



\## Requirements



\- ArcGIS Pro with ArcPy

\- Python environment associated with ArcGIS Pro

\- pandas

\- NumPy

\- Microsoft Excel support for reading the FVS species crosswalk

\- Forest Vegetation Simulator

\- TreeMap

\- FIA TREE and COND tables



\## Data availability



Large input datasets are not distributed with this repository. Users must

download TreeMap, FIA, ownership, and FVS reference data from their respective

official sources.



\## Status



LETO is under active development

