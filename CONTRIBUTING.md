# Contributing to LETO

## Development checks without ArcGIS

The package declares Python 3.9 or newer. These commands install the lightweight
package and test dependencies; they do not install ArcPy or the full GIS runtime.

Windows PowerShell, from the repository root:

```powershell
cd LETO_v1_project
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe -m pytest
```

macOS/Linux, for the unit tests:

```bash
cd LETO_v1_project
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pytest
```

Passing unit tests does not validate an end-to-end spatial run. Read
[AGENTS.md](AGENTS.md) before contributing; preserve reproducible regional configurations and document
any intentional changes to scientific assumptions.

When reporting a problem, include the failing stage, traceback, environment
versions, relevant configuration with private paths removed, and available QA
summaries. Do not attach large source datasets.
