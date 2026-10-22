"""Maine baseline settings with cluster paths; set LETO_DATA_ROOT and LETO_RUNS_ROOT."""
import importlib.util
import os
from pathlib import Path
from uuid import uuid4

spec = importlib.util.spec_from_file_location('maine_reference', Path(__file__).with_name('maine_test.py'))
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)
CONFIG = reference.CONFIG.copy()
data = Path(os.environ['LETO_DATA_ROOT'])
runs = Path(os.environ['LETO_RUNS_ROOT'])
job = os.environ.get('SLURM_JOB_ID', uuid4().hex)
CONFIG.update({
    'working_gdb': str(runs / f'working_{job}.workspace'),
    'runs_root_folder': str(runs),
    'study_area': str(data / 'Maine_5_co.gdb' / 'ME_5_Co_Dissolve'),
    'treemap_raster': str(data / 'TreeMap2022_CONUS.tif'),
    'fia_tree_files': [str(data / 'FIA_Tree' / f'{state}_TREE.csv') for state in ('MA', 'ME', 'NH', 'NY')],
    'fia_cond_files': [str(data / 'FIA_Cond' / f'{state}_COND.csv') for state in ('MA', 'ME', 'NH', 'NY')],
    'species_crosswalk': str(data / 'FVS_SpeciesCrosswalk.xls'),
    'blank_fvs_database': str(data / 'FVS_Data.db'),
    'ownership_source_raster': str(data / 'US_forest_ownership.tif'),
    'riparian': {
        **reference.CONFIG['riparian'],
        'nhd_flowline': str(data / 'NHDPlus_H_National_Release_2_GDB.gdb' / 'NHDPlus_H_National' / 'NetworkNHDFlowline'),
        'nhd_area': str(data / 'NHDPlus_H_National_Release_2_GDB.gdb' / 'NHDPlus_H_National' / 'NHDArea'),
        'nhd_waterbody': str(data / 'NHDPlus_H_National_Release_2_GDB.gdb' / 'NHDPlus_H_National' / 'NHDWaterbody'),
    },
})
