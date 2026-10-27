# Using LETO

This folder contains the application you run for your study area.
Start with the [repository overview](../README.md) for requirements and expected
outputs.

## Your first run

1. Follow [Getting started](docs/getting_started.md) to prepare ArcGIS Pro and
   your source datasets.
2. Copy `configs/maine_test.py` to `configs/my_region.py` and edit the copy.
   Use the [configuration guide](docs/configuration.md) to review your choices.
3. From this folder in an ArcGIS Pro Python Command Prompt, run:

```bat
python run_leto.py configs\my_region.py
```

For notebook execution, see the [setup guide](docs/getting_started.md).

## Use and review your outputs

- [Workflow and outputs](docs/workflow.md): locate management units, inventory
  tables, databases, and QA reports.
- [Validation](docs/validation.md): review the results before downstream use.
- [Troubleshooting](docs/troubleshooting.md): diagnose setup and data problems.
- [Advanced behavior](docs/advanced_behavior.md): understand small-stand
  handling, grouping, identifiers, and imputation diagnostics.

Use a dedicated working geodatabase and configuration for each independent
project. Save your configuration, data versions, and spatial outputs so you
can reproduce a result later.
