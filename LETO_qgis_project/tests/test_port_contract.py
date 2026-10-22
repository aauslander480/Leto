"""Guard the original scientific functions against accidental edits during porting."""
import ast
from pathlib import Path

from leto.gis_backend import _distance_meters, EnvManager, env


def test_distance_units_keep_foot_buffers():
    assert _distance_meters('45 Feet') == 13.716000000000001
    assert _distance_meters('60 Meters') == 60
    assert _distance_meters('2', 0.3048) == 0.6096


def test_scientific_function_bodies_are_unchanged():
    root = Path(__file__).parents[1]
    original = root.parent / 'LETO_v1_project' / 'leto'
    checked = []
    for path in sorted(original.glob('*.py')):
        old = ast.parse(path.read_text(encoding='utf-8'))
        new = ast.parse((root / 'leto' / path.name).read_text(encoding='utf-8'))
        functions = {n.name: n for n in new.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
        for node in old.body:
            if not isinstance(node, ast.FunctionDef):
                continue
            # Functions touching GIS are exercised separately with real PyQGIS.
            if any(isinstance(n, ast.Name) and n.id == 'arcpy' for n in ast.walk(node)):
                continue
            if node.name in {'run_pipeline', 'ensure_working_gdb'}:
                continue
            assert ast.dump(node) == ast.dump(functions[node.name]), (path.name, node.name)
            checked.append((path.name, node.name))
    assert len(checked) >= 70, f'Unexpectedly few preserved scientific functions: {len(checked)}'


def test_environment_restores_on_failure():
    old = env.snapRaster
    try:
        with EnvManager(snapRaster='temporary'):
            assert env.snapRaster == 'temporary'
            raise ValueError('test')
    except ValueError:
        pass
    assert env.snapRaster == old
