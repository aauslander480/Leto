# AGENTS.md

## Purpose

This repository contains LETO, a forest landscape segmentation and FVS inventory-generation workflow.

Codex and other coding agents must treat the verified Maine baseline as the source of truth unless a task explicitly requests a scientific change.

## Verified Baseline

The current verified baseline produces:

- Valid segmentation cells: 429,977
- Parent segments: 12,168
- Final management units: 19,662
- Direct FVS stands: 19,244
- Nearest-imputed stands: 418
- Grouped FVS databases: 6

When a task is not intended to change scientific behavior, these values are acceptance criteria.

## General Development Rules

1. Read this file, `README.md`, and `docs/baseline_validation.md` before editing code.
2. Preserve scientific behavior unless the task explicitly requests a scientific change.
3. Make one substantive feature change at a time.
4. Keep the verified baseline reproducible through configuration.
5. Never remove or overwrite a known-good workflow without preserving a rollback path.
6. Prefer small, reviewable commits with clear messages.
7. Do not silently alter thresholds, weights, classifications, state codes, species mappings, buffer distances, donor logic, or FVS sampling assumptions.
8. Do not introduce fixed local paths into package code.
9. Keep run-specific paths and modeling assumptions in regional configuration files.
10. Preserve backward compatibility when practical, especially for existing public entry points and configuration keys.

## Configuration Rules

All user-selectable modeling assumptions must be represented in configuration rather than hard-coded in stage modules.

This includes:

- Study-area and source-data paths
- TreeMap inputs
- FIA TREE and COND inputs
- Ownership settings
- Segmentation thresholds
- Objective-function weights
- Neighborhood and convergence settings
- Small-stand handling
- Similar-stand merge criteria
- Proximity-merge criteria
- Riparian enablement
- NHD flowline, area, and waterbody buffer rules
- Wetland settings
- Donor filtering
- Species crosswalk settings
- Missing-stand imputation
- FVS variant, state code, inventory year, and sampling design
- Database grouping behavior

Every resolved setting must be written to `Run_Settings.json`.

Configuration validation must fail early with clear messages when required values, paths, fields, or rule structures are missing or invalid.

## Scientific Integrity Rules

### Acreage

- Prefer preserving acreage over deleting stands.
- Any operation that removes cells or polygons must be explicit, configurable, documented, and reported.
- Report acreage before and after segmentation, merging, riparian splitting, and final management-unit creation.
- Unexpected acreage loss must raise an error or a clearly visible warning.

### Identifiers

Use identifiers consistently:

- `SEG_ID`: original physically connected biological segment
- `MU_ID`: final management unit; may represent one or more `SEG_ID` values
- `ParentSEG_ID` or `PARENT_SEG`: original segment provenance
- `MU_PART`: optional component number for multipart management units
- `REGION_ID`: regional identifier used to avoid cross-state collisions
- Global FVS stand identifiers should include the region when outputs may later be combined

Do not repurpose an existing identifier without updating documentation, crosswalk logic, QA, and tests.

### Donor Weights

- Donor weights must sum to 1 within each final `MU_ID` within numerical tolerance.
- If multiple segments are grouped into one management unit, combine donor weights using contributing segment acreage.
- Preserve donor provenance.
- Filtering donors must retain at least one donor per stand.
- Re-normalize retained weights after filtering.

### FIA Trees and Species

- Keep only intended FIA tree records, normally live trees unless configuration says otherwise.
- Preserve `PLT_CN` as a non-floating identifier.
- Do not silently drop untranslated species without QA output.
- Report record counts, weighted TPA, and affected stands for untranslated species.
- Preserve height fallback behavior unless explicitly changed.

### FVS Inventory

- Every tree `STAND_ID` must exist in the stand table.
- Every runnable stand must contain at least one valid tree record.
- `TREE_COUNT` must remain positive and represent the intended weighted TPA.
- Tree IDs must be unique within each stand.
- Database schema and sampling-design values must match the configured blank FVS database and documented workflow.

## Pipeline Structure

LETO currently follows five main stages:

1. Build segmentation attribute rasters
2. Segment the landscape and create management units
3. Build the segment/management-unit donor crosswalk
4. Build weighted FVS stand and tree inventories
5. Create FVS-ready SQLite databases

Keep stage responsibilities clear.

Prefer DataFrame handoffs between Stages 3, 4, and 5. CSV files should remain available as audit outputs but should not be required as internal handoffs when in-memory objects are available.

## Module Design

- Keep modules focused on one responsibility.
- Avoid large functions that combine reading, transformation, QA, formatting, and writing.
- Move shared validation, logging, identifier cleaning, and raster utility functions into shared modules.
- Public functions should have docstrings describing inputs, outputs, assumptions, and failure conditions.
- Avoid hidden dependencies on module-level mutable globals.
- Pass configuration or explicit parameters into functions.

## ArcPy and Workspace Rules

- Use isolated run folders and workspaces.
- Do not share temporary datasets between independent regional runs.
- Use `arcpy.EnvManager` where possible so environment changes are scoped.
- Use nearest-neighbor resampling for categorical rasters.
- Validate raster extent, dimensions, cell size, alignment, and coordinate system before cell-by-cell operations.
- Avoid simultaneous writes to the same geodatabase.
- Do not assume parallel ArcPy execution is safe until isolated-process testing is complete.

## Multi-Region Scalability

The software must remain capable of processing multiple states or regions independently.

Design requirements:

- One codebase may serve many regional configuration files.
- Each region must have its own `REGION_ID`, run folder, working datasets, logs, and outputs.
- A failure in one regional run must not invalidate completed runs.
- Avoid state-specific logic in core modules.
- State-specific buffer rules, FVS settings, and source paths belong in configuration or rule files.
- Sequential batch execution should be supported before parallel execution is attempted.

## Proximity-Merge Feature Rules

When proximity merging is implemented:

- It must be optional.
- Existing actions such as `force_merge`, `retain`, and `delete` must continue to work unchanged.
- The default baseline configuration must remain reproducible.
- Candidate units do not need to touch physically.
- Do not alter geometry to fill gaps.
- Preserve original `SEG_ID` values.
- Assign a shared `MU_ID` to grouped components.
- Require or evaluate ownership, forest type, biological similarity, distance, and combined acreage according to configuration.
- Retain unmatched small units when configured to do so.
- Report number of merges, distances, multipart-unit counts, remaining small units, acreage, and donor-weight validation.

## Testing Requirements

At minimum, maintain tests for:

- Configuration loading and validation
- Riparian-rule validation
- Raster alignment validation
- Identifier normalization
- Donor weights summing to 1
- Area-weighted aggregation across multiple parent segments
- Stand/tree referential integrity
- Unique tree IDs within stands
- FVS database schema and row counts
- Baseline mode preserving verified Maine outputs where practical

Scientific feature changes must include focused tests for the new behavior.

## QA and Reporting

Every run should record:

- LETO version
- Region ID
- Full resolved configuration
- Input paths
- Stage timings
- Cell, segment, management-unit, stand, tree, and donor counts
- Acreage totals
- Ownership and management-class summaries
- Missing-species diagnostics
- Missing or imputed stand diagnostics
- Database group summaries
- Output paths
- Warnings and deviations from expected behavior

Do not remove existing QA outputs without an explicit replacement.

## Change Management

Before editing:

1. Inspect the relevant modules and tests.
2. State which scientific behavior, if any, will change.
3. Identify baseline outputs that should remain unchanged.

After editing:

1. Run focused unit tests.
2. Run import and syntax checks.
3. Run the Maine baseline when ArcPy and source data are available.
4. Compare outputs with the verified baseline.
5. Summarize changed files, behavior, tests, and any unresolved risks.

## Commit Guidance

Use descriptive commits, for example:

- `Move riparian buffer rules into regional configuration`
- `Add validation for donor weight normalization`
- `Implement optional 60 m proximity grouping`
- `Preserve baseline small-stand handling`

Do not combine unrelated refactors and scientific changes in one commit.

## Current Priority

The current development priority is:

1. Preserve the verified baseline.
2. Keep all modeling assumptions in configuration.
3. Add the optional proximity-merge feature as an isolated, testable change.
4. Maintain independent multi-region scalability.
