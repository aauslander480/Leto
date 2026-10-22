"""Real GIS integration tests; run with QGIS's Python, never a mocked provider."""
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('qgis.core')
pytest.importorskip('osgeo.gdal')

from leto import gis_backend as gis


@pytest.fixture(autouse=True)
def runtime():
    gis.initialize()
    with gis.EnvManager(outputCoordinateSystem=None, snapRaster=None):
        yield
    gis._aliases.clear()


def vector(path, rectangles):
    q = gis.q
    layer = q.QgsVectorLayer('Polygon?crs=EPSG:5070&field=MU_ID:integer&field=SEG_ID:integer', 'fixture', 'memory')
    for ident, bounds in rectangles:
        feature = q.QgsFeature(layer.fields())
        feature.setGeometry(q.QgsGeometry.fromRect(q.QgsRectangle(*bounds)))
        feature.setAttributes([ident, ident])
        assert layer.dataProvider().addFeatures([feature])[0]
    gis._save(layer, str(path))
    return str(path)


def raster(path, array, nodata=-9999):
    gis.NumPyArrayToRaster(array, gis.Point(0, 0), 30, 30, nodata).save(str(path))
    crs = gis.q.QgsCoordinateReferenceSystem('EPSG:5070')
    gis.management.DefineProjection(str(path), gis.SpatialReference(crs))
    return str(path)


def test_raster_orientation_nodata_rat_and_polygonization(tmp_path):
    array = np.array([[1, -9999, 2], [1, 1, 2]], dtype=np.int32)
    source = raster(tmp_path / 'source.tif', array)
    assert np.array_equal(gis.RasterToNumPyArray(source, nodata_to_value=-9), [[1, -9, 2], [1, 1, 2]])
    desc = gis.Describe(source)
    assert (desc.extent.YMin, desc.extent.YMax, desc.width, desc.height) == (0, 60, 3, 2)
    gis.management.BuildRasterAttributeTable(source)
    with gis.da.SearchCursor(source, ['VALUE', 'COUNT']) as cursor:
        assert list(cursor) == [(1, 3), (2, 2)]
    polygons = str(tmp_path / 'polygons')
    gis.conversion.RasterToPolygon(source, polygons, 'NO_SIMPLIFY', 'Value', 'SINGLE_OUTER_PART')
    assert gis.management.GetCount(polygons)[0] == 2
    layer = gis._vector(polygons)
    assert sum(f.geometry().area() for f in layer.getFeatures()) == 5 * 900


def test_clipping_preserves_source_grid_and_center_mask(tmp_path):
    source = raster(tmp_path / 'source.tif', np.arange(9, dtype=np.int32).reshape(3, 3))
    boundary = vector(tmp_path / 'boundary', [(1, (5, 5, 55, 55))])
    clipped = str(tmp_path / 'clipped')
    gis.management.Clip(source, '#', clipped, boundary)
    desc = gis.Describe(clipped)
    assert (desc.extent.XMin, desc.extent.YMin, desc.width, desc.height) == (0, 0, 2, 2)
    assert np.array_equal(gis.RasterToNumPyArray(clipped), [[3, 4], [6, 7]])


def test_categorical_alignment_masks_reference_nodata(tmp_path):
    source = raster(tmp_path / 'ownership.tif', np.array([[3, 4], [5, 6]], dtype=np.int16))
    reference = raster(tmp_path / 'reference.tif', np.array([[1, -9999], [1, 1]], dtype=np.int32))
    output = str(tmp_path / 'aligned')
    gis.sa.ExtractByMask(source, reference).save(output)
    assert np.array_equal(gis.RasterToNumPyArray(output, nodata_to_value=0), [[3, 0], [5, 6]])


def test_cursors_join_delete_and_geodesic_area(tmp_path):
    path = vector(tmp_path / 'units', [(1, (0, 0, 30, 30)), (2, (60, 0, 90, 30))])
    gis.management.AddField(path, 'Acres', 'DOUBLE')
    with gis.da.UpdateCursor(path, ['MU_ID', 'Acres', 'SHAPE@']) as cursor:
        for row in cursor:
            if row[0] == 2:
                cursor.deleteRow()
            else:
                row[1] = row[2].getArea('GEODESIC', 'ACRES')
                cursor.updateRow(row)
    with gis.da.SearchCursor(path, ['MU_ID', 'Acres']) as cursor:
        rows = list(cursor)
    assert len(rows) == 1 and rows[0][0] == 1 and 0.2 < rows[0][1] < 0.25
    table = tmp_path / 'summary.csv'
    table.write_text('SEG_ID,AGE\n1,22\n')
    gis.management.JoinField(path, 'SEG_ID', str(table), 'SEG_ID', ['AGE'])
    with gis.da.SearchCursor(path, ['AGE']) as cursor:
        assert list(cursor) == [[22]]


def test_near_uses_polygon_boundaries_and_preserves_fids(tmp_path):
    path = vector(tmp_path / 'units', [(1, (0, 0, 10, 10)), (2, (20, 0, 30, 10)), (3, (50, 0, 60, 10))])
    gis.management.MakeFeatureLayer(path, 'missing')
    gis.management.MakeFeatureLayer(path, 'donors')
    gis.management.SelectLayerByAttribute('missing', 'NEW_SELECTION', 'MU_ID = 1')
    gis.management.SelectLayerByAttribute('donors', 'NEW_SELECTION', 'MU_ID IN (2, 3)')
    output = str(tmp_path / 'near')
    gis.analysis.GenerateNearTable('missing', 'donors', output, closest_count=1)
    with gis.da.SearchCursor(output, ['IN_FID', 'NEAR_FID', 'NEAR_DIST']) as cursor:
        rows = list(cursor)
    with gis.da.SearchCursor(path, ['OID@', 'MU_ID']) as cursor:
        ids = {mu: oid for oid, mu in cursor}
    assert rows == [(ids[1], ids[2], 10)]


def test_near_ties_stop_and_export_candidates(tmp_path):
    path = vector(tmp_path / 'units', [(1, (0, 0, 10, 10)), (2, (20, 0, 30, 10)), (3, (-20, 0, -10, 10))])
    gis.management.MakeFeatureLayer(path, 'missing')
    gis.management.MakeFeatureLayer(path, 'donors')
    gis.management.SelectLayerByAttribute('missing', 'NEW_SELECTION', 'MU_ID = 1')
    gis.management.SelectLayerByAttribute('donors', 'NEW_SELECTION', 'MU_ID IN (2, 3)')
    with pytest.raises(RuntimeError, match='ties'):
        gis.analysis.GenerateNearTable('missing', 'donors', str(tmp_path / 'near'))
    assert (tmp_path / 'near.ties.csv').exists()


def test_overlay_buffers_keep_acreage(tmp_path):
    source = vector(tmp_path / 'source', [(1, (0, 0, 100, 100))])
    overlay = vector(tmp_path / 'overlay', [(2, (40, 0, 60, 100))])
    intersect, erase = str(tmp_path / 'intersect'), str(tmp_path / 'erase')
    gis.analysis.PairwiseIntersect([source, overlay], intersect)
    gis.analysis.PairwiseErase(source, overlay, erase)
    area = lambda p: sum(f.geometry().area() for f in gis._vector(p).getFeatures())
    assert area(intersect) + area(erase) == 10000
    buffered = str(tmp_path / 'buffer')
    gis.analysis.PairwiseBuffer(overlay, buffered, '45 Feet', dissolve_option='ALL')
    bounds = gis._vector(buffered).extent()
    assert bounds.xMinimum() == pytest.approx(40 - 13.716)
    dissolved = str(tmp_path / 'dissolved')
    gis.management.Dissolve(source, dissolved, 'MU_ID')
    assert area(dissolved) == 10000


def test_geodesic_near_is_never_silently_planar(tmp_path):
    with pytest.raises(NotImplementedError, match='GEODESIC'):
        gis.analysis.GenerateNearTable('unused', 'unused', str(tmp_path / 'near'), method='GEODESIC')


def test_stage_polygon_finalization_and_donor_handoff(tmp_path, monkeypatch):
    """Exercise real stage calls together, beyond individual backend methods."""
    import pandas as pd
    from leto import stage2_segmentation as s2, stage3_crosswalk as s3
    from leto.pipeline import _ensure_management_units
    from leto.spatial_management_units import finalize_management_outputs
    labels = raster(tmp_path / 'labels.tif', np.array([[1, 2]], dtype=np.int32))
    treemap = raster(tmp_path / 'treemap.tif', np.array([[10, 20]], dtype=np.int32))
    polygons = str(tmp_path / 'segments')
    summary = pd.DataFrame({'SEG_ID': [1, 2], 'Acres': [900/4046.8564224]*2,
                            'STDAGE_MEAN': [30., 40.], 'FORTYPCD_DOM': [101, 101],
                            'OWN_TYPE': ['Family', 'Family'], 'OWN_CODE': [3, 3]})
    csv = tmp_path / 'summary.csv'
    summary.to_csv(csv, index=False)
    monkeypatch.setattr(s2, 'segment_raster', labels)
    monkeypatch.setattr(s2, 'segment_polygons', polygons)
    monkeypatch.setattr(s2, 'segment_summary_csv', str(csv))
    monkeypatch.setattr(s2, 'output_gdb', str(tmp_path))
    s2.write_segment_polygons(summary)
    components, units = str(tmp_path / 'components'), str(tmp_path / 'units')
    _ensure_management_units(polygons, components)
    qa = finalize_management_outputs(components, components, units,
         proximity_enabled=False, maximum_distance_meters=60,
         warning_chain_ratio=2, qa_folder=str(tmp_path), acreage_tolerance_acres=0.01)
    assert qa['component_polygons'] == qa['final_management_units'] == 2
    info = s3.validate_raster_alignment(labels, treemap)
    parents = s3.build_weights_table(labels, treemap, {10: 'A', 20: 'B'}, s3.cell_area_acres(info))
    stands = s3.read_management_units(units)
    parts = s3.read_management_segments(components)
    weights = s3.expand_parent_weights_to_management_units(parents, stands, parts)
    assert set(weights['PLT_CN']) == {'A', 'B'}
    assert all(abs(v - 1) < 1e-12 for v in weights.groupby('MU_ID')['WEIGHT'].sum())
