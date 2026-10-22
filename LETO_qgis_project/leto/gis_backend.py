"""LETO's spatial I/O contract implemented with PyQGIS and GDAL.

The small, explicit API retains the original stage call signatures so the
scientific array and table algorithms need no rewrite. This is not an ArcPy
emulator: only operations used by LETO are implemented. Unsupported semantics
raise errors. Outputs are GeoTIFFs and individual GeoPackages in a workspace
directory; source file geodatabases remain read-only.
"""
from __future__ import annotations

import math
import os
import re
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

env = SimpleNamespace(overwriteOutput=True, snapRaster=None, cellSize=None,
                      extent=None, outputCoordinateSystem=None)
ExecuteError = RuntimeError
_aliases = {}
_app = None
_ready = False
BUFFER_SEGMENTS = 64  # segments per quadrant; recorded in Backend_Settings.json


def initialize():
    """Initialize a headless QGIS session using its installed Python runtime."""
    global _app, _ready, q, gdal, ogr, osr, QVariant, processing
    if _ready:
        return
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    try:
        from qgis import core as q
        from qgis.PyQt.QtCore import QVariant
        from osgeo import gdal, ogr, osr
    except ImportError as exc:
        raise RuntimeError('Use the Python supplied by the QGIS module; PyQGIS '
                           'and GDAL bindings are required. Ordinary pip Python is insufficient.') from exc
    if q.QgsApplication.instance() is None:
        if os.environ.get('QGIS_PREFIX_PATH'):
            q.QgsApplication.setPrefixPath(os.environ['QGIS_PREFIX_PATH'], True)
        _app = q.QgsApplication([], False)
        _app.initQgis()
    import sys
    plugins = str(Path(q.QgsApplication.pkgDataPath()) / 'python' / 'plugins')
    if plugins not in sys.path:
        sys.path.insert(0, plugins)
    import processing
    from processing.core.Processing import Processing
    Processing.initialize()
    gdal.UseExceptions()
    _ready = True


def AddMessage(text):
    # Stage functions already print their own messages.
    pass


@contextmanager
def EnvManager(**kwargs):
    old = {key: getattr(env, key, None) for key in kwargs}
    for key, value in kwargs.items():
        setattr(env, key, value)
    try:
        yield
    finally:
        for key, value in old.items():
            setattr(env, key, value)


def _file(path, kind=None):
    path = str(path)
    if Path(path).suffix.lower() in {'.tif', '.tiff', '.gpkg', '.csv', '.shp', '.dbf'}:
        return path
    if kind:
        return path + {'raster': '.tif', 'vector': '.gpkg', 'table': '.csv'}[kind]
    for suffix in ('', '.tif', '.gpkg', '.csv'):
        if Path(path + suffix).exists():
            return path + suffix
    return path


def _uri(path):
    path = _file(path)
    # ArcGIS feature-dataset components are not part of the OGR layer name.
    match = re.match(r'^(.*?\.gdb)[/\\](.+)$', path, re.I)
    if match:
        return match[1] + '|layername=' + re.split(r'[/\\]', match[2])[-1]
    return path


def _vector(path):
    initialize()
    if isinstance(path, q.QgsVectorLayer):
        return path
    if str(path) in _aliases:
        return _aliases[str(path)]
    layer = q.QgsVectorLayer(_uri(path), Path(str(path)).stem, 'ogr')
    if not layer.isValid():
        raise RuntimeError(f'Cannot open vector dataset: {path}')
    return layer


def _raster(path, update=False):
    initialize()
    ds = gdal.Open(_file(path), gdal.GA_Update if update else gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f'Cannot open raster: {path}')
    if ds.RasterCount != 1:
        raise ValueError(f'LETO requires a single-band raster: {path}')
    return ds


def Exists(path):
    if str(path) in _aliases or Path(_file(path)).exists():
        return True
    if '.gdb' in str(path).lower() or '|layername=' in str(path):
        try:
            _vector(path)
            return True
        except RuntimeError:
            pass
    return False


class SpatialReference:
    def __init__(self, crs):
        initialize()
        self.crs = crs
        if not crs.isValid():
            raise ValueError('Dataset has no valid coordinate reference system.')
        self.type = 'Geographic' if crs.isGeographic() else 'Projected'
        self.factoryCode = int(crs.authid().split(':')[-1]) if crs.authid().startswith('EPSG:') else 0
        self.linearUnitName = q.QgsUnitTypes.toString(crs.mapUnits())
        self.metersPerUnit = (q.QgsUnitTypes.fromUnitToUnitFactor(
            crs.mapUnits(), q.QgsUnitTypes.DistanceMeters) if not crs.isGeographic() else None)

    def exportToString(self):
        return self.crs.toWkt()


class Point:
    def __init__(self, x, y):
        self.X, self.Y = float(x), float(y)


class Geometry:
    def __init__(self, geometry, crs):
        self.geometry, self.crs = geometry, crs

    @property
    def trueCentroid(self):
        point = self.geometry.centroid().asPoint()
        return Point(point.x(), point.y())

    def distanceTo(self, other):
        if self.crs != other.crs:
            raise ValueError('Distance requires matching coordinate systems.')
        return self.geometry.distance(other.geometry)

    def getArea(self, method, units):
        if method != 'GEODESIC' or units != 'ACRES':
            raise NotImplementedError((method, units))
        measure = q.QgsDistanceArea()
        measure.setSourceCrs(self.crs, q.QgsProject.instance().transformContext())
        ellipsoid = self.crs.ellipsoidAcronym()
        if not ellipsoid or not measure.setEllipsoid(ellipsoid):
            raise ValueError('Cannot determine the dataset ellipsoid for geodesic acreage.')
        return measure.measureArea(self.geometry) / 4046.8564224


def PointGeometry(point, spatial_reference):
    initialize()
    return Geometry(q.QgsGeometry.fromPointXY(q.QgsPointXY(point.X, point.Y)), spatial_reference.crs)


def Describe(path):
    initialize()
    if Path(_file(path)).suffix.lower() in {'.tif', '.tiff'}:
        ds = _raster(path)
        gt = ds.GetGeoTransform()
        if gt[2] != 0 or gt[4] != 0 or gt[1] <= 0 or gt[5] >= 0:
            raise ValueError('LETO requires a north-up, unrotated raster.')
        crs = q.QgsCoordinateReferenceSystem.fromWkt(ds.GetProjection())
        extent = SimpleNamespace(XMin=gt[0], YMin=gt[3] + ds.RasterYSize * gt[5],
                                 XMax=gt[0] + ds.RasterXSize * gt[1], YMax=gt[3])
        return SimpleNamespace(width=ds.RasterXSize, height=ds.RasterYSize,
                               meanCellWidth=gt[1], meanCellHeight=-gt[5], extent=extent,
                               spatialReference=SpatialReference(crs))
    layer = _vector(path)
    return SimpleNamespace(spatialReference=SpatialReference(layer.crs()))


def Raster(path):
    return Describe(path)


def RasterToNumPyArray(path, nodata_to_value=None):
    ds = _raster(path)
    band = ds.GetRasterBand(1)
    array = band.ReadAsArray()
    if nodata_to_value is not None:
        valid = band.GetMaskBand().ReadAsArray() != 0
        dtype = np.result_type(array.dtype, np.min_scalar_type(nodata_to_value))
        array = array.astype(dtype)
        array[~valid] = nodata_to_value
    return array


def _write_raster(path, array, gt, projection, nodata):
    initialize()
    from osgeo import gdal_array
    dest = _file(path, 'raster')
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    if array.dtype == np.dtype('bool'):
        array = array.astype(np.uint8)
    dtype = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype)
    if not dtype:
        raise TypeError(f'Unsupported raster dtype: {array.dtype}')
    ds = gdal.GetDriverByName('GTiff').Create(dest, array.shape[1], array.shape[0], 1, dtype,
                                            options=['COMPRESS=LZW', 'BIGTIFF=IF_SAFER'])
    ds.SetGeoTransform(gt)
    if projection:
        ds.SetProjection(projection)
    band = ds.GetRasterBand(1)
    if nodata is not None:
        band.SetNoDataValue(float(nodata))
    band.WriteArray(array)
    ds.FlushCache()
    ds = None
    return dest


class _ArrayRaster:
    def __init__(self, array, point, dx, dy, nodata, projection=''):
        self.array = array
        self.gt = (point.X, float(dx), 0, point.Y + array.shape[0] * float(dy), 0, -float(dy))
        self.nodata, self.projection = nodata, projection

    def save(self, path):
        return _write_raster(path, self.array, self.gt, self.projection, self.nodata)


def NumPyArrayToRaster(array, lower_left, x_cell_size, y_cell_size, value_to_nodata=None):
    return _ArrayRaster(array, lower_left, x_cell_size, y_cell_size, value_to_nodata)


def _rat_frame(path):
    """Read the original TreeMap RAT, including its PLT_CN and fallback columns."""
    ds = _raster(path)
    rat = ds.GetRasterBand(1).GetDefaultRAT()
    if rat is not None and rat.GetRowCount():
        columns = {}
        for col in range(rat.GetColumnCount()):
            reader = (rat.GetValueAsString if rat.GetTypeOfCol(col) == gdal.GFT_String
                      else rat.GetValueAsInt if rat.GetTypeOfCol(col) == gdal.GFT_Integer
                      else rat.GetValueAsDouble)
            columns[rat.GetNameOfCol(col)] = [reader(row, col) for row in range(rat.GetRowCount())]
        return pd.DataFrame(columns)
    vat = str(_file(path)) + '.vat.dbf'
    if Path(vat).exists():
        layer = _vector(vat)
        names = [field.name() for field in layer.fields()]
        return pd.DataFrame([[_clean(f[n]) for n in names] for f in layer.getFeatures()], columns=names)
    raise RuntimeError(f'Raster attribute table missing for {path}. Copy the original '
                       '.tif.vat.dbf and .aux.xml sidecars with TreeMap; a VALUE-only RAT is insufficient.')


def _clean(value):
    if value is None or (hasattr(value, 'isNull') and value.isNull()):
        return None
    return value


def ListFields(path):
    initialize()
    suffix = Path(_file(path)).suffix.lower()
    if suffix in {'.tif', '.tiff', '.csv'}:
        frame = _rat_frame(path) if suffix != '.csv' else pd.read_csv(_file(path))
        return [SimpleNamespace(name=str(n), type=('Double' if pd.api.types.is_numeric_dtype(frame[n]) else 'String'),
                                length=255) for n in frame.columns]
    layer = _vector(path)
    kinds = {QVariant.Int: 'Integer', QVariant.LongLong: 'BigInteger', QVariant.Double: 'Double',
             QVariant.String: 'String', QVariant.Date: 'Date'}
    return [SimpleNamespace(name=f.name(), type=('OID' if f.name().lower() == 'fid' else kinds.get(f.type(), 'String')),
                            length=f.length()) for f in layer.fields()]


class _Cursor:
    def __init__(self, path, fields, update=False):
        self.path, self.fields, self.update = path, list(fields), update
        self.layer = None

    def __enter__(self):
        suffix = Path(_file(self.path)).suffix.lower()
        if suffix in {'.tif', '.tiff', '.csv'}:
            if self.update:
                raise NotImplementedError('Updating raster RATs or CSVs through cursors is not supported.')
            frame = _rat_frame(self.path) if suffix != '.csv' else pd.read_csv(_file(self.path))
            self.iterator = iter(frame[self.fields].itertuples(index=False, name=None))
        else:
            self.layer = _vector(self.path)
            missing = [n for n in self.fields if n not in {'OID@', 'SHAPE@'} and self.layer.fields().indexFromName(n) < 0]
            if missing:
                raise ValueError(f'Missing fields in {self.path}: {missing}')
            if self.update and not self.layer.startEditing():
                raise RuntimeError(f'Cannot edit {self.path}')
            self.iterator = iter(self.layer.getSelectedFeatures() if str(self.path) in _aliases
                                 else self.layer.getFeatures())
        return self

    def __iter__(self):
        return self

    def __next__(self):
        row = next(self.iterator)
        if self.layer is None:
            return row
        self.current = row
        return [row.id() if n == 'OID@' else Geometry(row.geometry(), self.layer.crs()) if n == 'SHAPE@'
                else _clean(row[n]) for n in self.fields]

    def updateRow(self, values):
        if not self.update:
            raise RuntimeError('Read-only cursor')
        for name, value in zip(self.fields, values):
            if name == 'OID@':
                continue
            if name == 'SHAPE@':
                if not self.layer.changeGeometry(self.current.id(), value.geometry):
                    raise RuntimeError('Geometry update failed')
            elif not self.layer.changeAttributeValue(self.current.id(), self.layer.fields().indexFromName(name), value):
                raise RuntimeError(f'Attribute update failed: {name}')

    def deleteRow(self):
        if not self.update or not self.layer.deleteFeature(self.current.id()):
            raise RuntimeError('Feature deletion failed')

    def __exit__(self, typ, value, traceback):
        if self.update:
            if typ:
                self.layer.rollBack()
            elif not self.layer.commitChanges():
                raise RuntimeError(f'Failed to commit {self.path}: {self.layer.commitErrors()}')
        self.iterator = None
        self.layer = None


da = SimpleNamespace(SearchCursor=lambda path, fields: _Cursor(path, fields),
                     UpdateCursor=lambda path, fields: _Cursor(path, fields, True))


def _save(layer, path):
    initialize()
    dest = _file(path, 'vector')
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    if Path(dest).exists():
        management.Delete(dest)
    options = q.QgsVectorFileWriter.SaveVectorOptions()
    options.driverName = 'GPKG'
    options.layerName = Path(dest).stem
    # Do not copy provider FIDs as ordinary attributes into a new GeoPackage.
    options.attributes = [i for i, f in enumerate(layer.fields()) if f.name().lower() != 'fid']
    result = q.QgsVectorFileWriter.writeAsVectorFormatV3(
        layer, dest, q.QgsProject.instance().transformContext(), options)
    if result[0] != q.QgsVectorFileWriter.NoError:
        raise RuntimeError(f'Writing {dest} failed: {result}')
    return dest


def _run(algorithm, **params):
    initialize()
    context = q.QgsProcessingContext()
    context.setInvalidGeometryCheck(q.QgsFeatureRequest.GeometryAbortOnInvalid)
    return processing.run(algorithm, {**params, 'OUTPUT': 'memory:'}, context=context)['OUTPUT']


def _model_crs():
    value = env.outputCoordinateSystem
    if value is None:
        return None
    if isinstance(value, SpatialReference):
        return value.crs
    return Describe(value).spatialReference.crs


def _projected(path, target=None):
    layer = _vector(path)
    target = target or _model_crs() or layer.crs()
    if not target.isValid() or target.isGeographic():
        raise ValueError('LETO vector processing requires the projected model CRS.')
    if not layer.crs().isValid():
        raise ValueError(f'Unknown input CRS: {path}')
    if layer.crs() != target:
        layer = _run('native:reprojectlayer', INPUT=layer, TARGET_CRS=target)
    return layer


class _Management:
    def create_workspace(self, parent, name):
        path = Path(parent) / name
        path.mkdir(parents=True, exist_ok=True)
        return str(path)

    def Delete(self, path):
        if str(path) in _aliases:
            del _aliases[str(path)]
            return
        dest = Path(_file(path))
        if dest.is_dir():
            raise ValueError('Dataset deletion cannot recursively delete a workspace.')
        if dest.exists():
            dest.unlink()
        for suffix in ('.aux.xml', '.ovr', '.msk', '-wal', '-shm'):
            side = Path(str(dest) + suffix)
            if side.exists():
                side.unlink()

    def CopyFeatures(self, source, output):
        return _save(_projected(source), output)

    def AddField(self, path, name, kind, field_length=255, **kwargs):
        layer = _vector(path)
        kinds = {'SHORT': QVariant.Int, 'LONG': QVariant.Int, 'BIGINTEGER': QVariant.LongLong,
                 'DOUBLE': QVariant.Double, 'FLOAT': QVariant.Double, 'TEXT': QVariant.String, 'DATE': QVariant.Date}
        if kind.upper() not in kinds:
            raise NotImplementedError(f'Field type: {kind}')
        field = q.QgsField(name, kinds[kind.upper()], len=field_length if kind.upper() == 'TEXT' else 0)
        if not layer.dataProvider().addAttributes([field]):
            raise RuntimeError(f'Cannot add field {name} to {path}')
        layer.updateFields()

    def CalculateField(self, path, field, expression, expression_type='PYTHON3'):
        # LETO uses field copies or Python literals only; no general eval.
        import ast
        match = re.fullmatch(r'!([^!]+)!', expression)
        source = match[1] if match else None
        literal = None if source else ast.literal_eval(expression)
        fields = [field, source] if source else [field]
        with da.UpdateCursor(path, fields) as cursor:
            for row in cursor:
                row[0] = row[1] if source else literal
                cursor.updateRow(row)

    def GetCount(self, path):
        layer = _vector(path)
        return [layer.selectedFeatureCount() if str(path) in _aliases else layer.featureCount()]

    def MakeFeatureLayer(self, source, name):
        # A separate provider instance keeps selections independent and FIDs stable.
        layer = _vector(source)
        _aliases[name] = q.QgsVectorLayer(layer.source(), name, layer.providerType())
        return name

    def SelectLayerByAttribute(self, name, selection_type, where_clause):
        if selection_type != 'NEW_SELECTION':
            raise NotImplementedError(selection_type)
        layer = _aliases[name]
        expression = q.QgsExpression(where_clause)
        if expression.hasParserError():
            raise ValueError(expression.parserErrorString())
        layer.selectByExpression(where_clause)

    def JoinField(self, path, key, table, table_key, fields):
        frame = pd.read_csv(_file(table))
        if frame[table_key].duplicated().any():
            raise ValueError('Join table keys must be unique.')
        lookup = frame.set_index(table_key).to_dict('index')
        for name in fields:
            if name not in [f.name for f in ListFields(path)]:
                dtype = frame[name].dtype
                self.AddField(path, name, 'LONG' if pd.api.types.is_integer_dtype(dtype)
                              else 'DOUBLE' if pd.api.types.is_numeric_dtype(dtype) else 'TEXT')
        with da.UpdateCursor(path, [key, *fields]) as cursor:
            for row in cursor:
                data = lookup.get(row[0], {})
                row[1:] = [None if pd.isna(data.get(n)) else data[n] for n in fields]
                cursor.updateRow(row)

    def Merge(self, sources, output):
        layers = [_projected(p) for p in sources]
        return _save(_run('native:mergevectorlayers', LAYERS=layers, CRS=layers[0].crs()), output)

    def MultipartToSinglepart(self, source, output):
        return _save(_run('native:multiparttosingleparts', INPUT=_projected(source)), output)

    def Dissolve(self, source, output, dissolve_field=None, multi_part='MULTI_PART'):
        if multi_part != 'MULTI_PART':
            raise NotImplementedError(multi_part)
        fields = [dissolve_field] if isinstance(dissolve_field, str) else dissolve_field or []
        return _save(_run('native:dissolve', INPUT=_projected(source), FIELD=fields), output)

    def DefineProjection(self, path, spatial_reference):
        ds = _raster(path, True)
        ds.SetProjection(spatial_reference.exportToString())
        ds.FlushCache()

    def CalculateStatistics(self, path):
        ds = _raster(path, True)
        ds.GetRasterBand(1).ComputeStatistics(False)
        ds.FlushCache()

    def BuildPyramids(self, path):
        ds = _raster(path, True)
        ds.BuildOverviews('NEAREST', [2, 4, 8, 16])
        ds.FlushCache()

    def BuildRasterAttributeTable(self, path, mode='Overwrite'):
        ds = _raster(path, True)
        band = ds.GetRasterBand(1)
        array = band.ReadAsArray()
        values, counts = np.unique(array[band.GetMaskBand().ReadAsArray() != 0], return_counts=True)
        rat = gdal.RasterAttributeTable()
        rat.CreateColumn('VALUE', gdal.GFT_Integer, gdal.GFU_MinMax)
        rat.CreateColumn('COUNT', gdal.GFT_Real, gdal.GFU_PixelCount)
        rat.SetRowCount(len(values))
        for i, (value, count) in enumerate(zip(values, counts)):
            rat.SetValueAsInt(i, 0, int(value))
            rat.SetValueAsDouble(i, 1, float(count))
        band.SetDefaultRAT(rat)
        ds.FlushCache()

    def Clip(self, in_raster, rectangle, out_raster, in_template_dataset, nodata_value='#',
             clipping_geometry='ClippingGeometry', maintain_clipping_extent='NO_MAINTAIN_EXTENT'):
        if clipping_geometry != 'ClippingGeometry' or maintain_clipping_extent != 'NO_MAINTAIN_EXTENT':
            raise NotImplementedError('Only source-grid-aligned polygon clipping is supported.')
        source = _raster(in_raster)
        gt = source.GetGeoTransform()
        sr = Describe(in_raster).spatialReference
        boundary = _vector(in_template_dataset)
        if boundary.crs() != sr.crs:
            boundary = _run('native:reprojectlayer', INPUT=boundary, TARGET_CRS=sr.crs)
        extent = boundary.extent()
        # Snap outward to source cell edges, matching NO_MAINTAIN_EXTENT.
        x0 = max(0, math.floor((extent.xMinimum() - gt[0]) / gt[1]))
        x1 = min(source.RasterXSize, math.ceil((extent.xMaximum() - gt[0]) / gt[1]))
        y0 = max(0, math.floor((gt[3] - extent.yMaximum()) / -gt[5]))
        y1 = min(source.RasterYSize, math.ceil((gt[3] - extent.yMinimum()) / -gt[5]))
        if x1 <= x0 or y1 <= y0:
            raise ValueError('Study boundary does not overlap TreeMap.')
        window_gt = (gt[0] + x0 * gt[1], gt[1], 0, gt[3] + y0 * gt[5], 0, gt[5])
        data = source.GetRasterBand(1).ReadAsArray(x0, y0, x1-x0, y1-y0)
        source_valid = source.GetRasterBand(1).GetMaskBand().ReadAsArray(x0, y0, x1-x0, y1-y0) != 0
        mask = _polygon_mask(boundary, window_gt, data.shape, source.GetProjection())
        nodata = source.GetRasterBand(1).GetNoDataValue()
        if nodata is None:
            nodata = -999999999
        data = data.astype(np.result_type(data.dtype, np.min_scalar_type(nodata)))
        data[~mask | ~source_valid] = nodata
        return _write_raster(out_raster, data, window_gt, source.GetProjection(), nodata)

    def ProjectRaster(self, source, output, spatial_reference, resampling, cell_size):
        if resampling != 'NEAREST':
            raise ValueError('Categorical ownership must use nearest-neighbor resampling.')
        return _align_raster(source, output, env.snapRaster)


def _polygon_mask(layer, gt, shape, projection):
    """Cell-center polygon mask; ALL_TOUCHED is deliberately disabled."""
    initialize()
    vector = ogr.GetDriverByName('Memory').CreateDataSource('')
    sr = osr.SpatialReference()
    sr.ImportFromWkt(projection)
    polygons = vector.CreateLayer('mask', srs=sr, geom_type=ogr.wkbUnknown)
    for feature in layer.getFeatures():
        geometry = feature.geometry()
        if not geometry.isGeosValid():
            raise ValueError('Invalid clipping geometry; repair the source explicitly.')
        row = ogr.Feature(polygons.GetLayerDefn())
        row.SetGeometry(ogr.CreateGeometryFromWkb(bytes(geometry.asWkb())))
        polygons.CreateFeature(row)
    mask = gdal.GetDriverByName('MEM').Create('', shape[1], shape[0], 1, gdal.GDT_Byte)
    mask.SetGeoTransform(gt)
    mask.SetProjection(projection)
    gdal.RasterizeLayer(mask, [1], polygons, burn_values=[1])
    return mask.ReadAsArray().astype(bool)


def _align_raster(source, output, reference):
    initialize()
    if reference is None:
        raise ValueError('An explicit reference raster is required for categorical alignment.')
    ref = _raster(reference)
    gt = ref.GetGeoTransform()
    dest = _file(output, 'raster')
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    ds = gdal.Warp(dest, _file(source), format='GTiff', dstSRS=ref.GetProjection(),
                   outputBounds=(gt[0], gt[3]+ref.RasterYSize*gt[5], gt[0]+ref.RasterXSize*gt[1], gt[3]),
                   width=ref.RasterXSize, height=ref.RasterYSize, resampleAlg='near',
                   outputType=gdal.GDT_Int32, dstNodata=-999999999,
                   creationOptions=['COMPRESS=LZW', 'BIGTIFF=IF_SAFER'])
    if ds is None:
        raise RuntimeError('Ownership reprojection failed.')
    ds.FlushCache()
    ds = None
    return dest


class _SpatialAnalyst:
    def ExtractByMask(self, source, reference):
        import tempfile
        with tempfile.TemporaryDirectory(prefix='leto_align_') as folder:
            dest = _align_raster(source, str(Path(folder) / 'aligned.tif'), reference)
            ds = _raster(dest)
            data = ds.ReadAsArray()
            ref = _raster(reference)
            valid = ref.GetRasterBand(1).GetMaskBand().ReadAsArray() != 0
            data[~valid] = -999999999
            gt, projection = ds.GetGeoTransform(), ds.GetProjection()
            ds = None
        return _ArrayRaster(data, Point(gt[0], gt[3]+data.shape[0]*gt[5]), gt[1], -gt[5],
                            -999999999, projection)


class _Conversion:
    def TableToTable(self, source, folder, name):
        dest = _file(str(Path(folder) / name), 'table')
        pd.read_csv(source).to_csv(dest, index=False)
        return dest

    def RasterToPolygon(self, in_raster, out_polygon_features, simplify, raster_field,
                        create_multipart_features):
        if simplify != 'NO_SIMPLIFY' or create_multipart_features != 'SINGLE_OUTER_PART':
            raise NotImplementedError('Only unsimplified singlepart polygonization is supported.')
        source = _raster(in_raster)
        sr = osr.SpatialReference()
        sr.ImportFromWkt(source.GetProjection())
        dest = _file(out_polygon_features, 'vector')
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        if Path(dest).exists():
            management.Delete(dest)
        vector = ogr.GetDriverByName('GPKG').CreateDataSource(dest)
        layer = vector.CreateLayer(Path(dest).stem, srs=sr, geom_type=ogr.wkbPolygon)
        layer.CreateField(ogr.FieldDefn('gridcode', ogr.OFTInteger64))
        band = source.GetRasterBand(1)
        # GDAL's default is four-connectivity, with no smoothing/simplification.
        gdal.Polygonize(band, band.GetMaskBand(), layer, 0, [])
        vector = None
        return dest


class _Analysis:
    def PairwiseDissolve(self, source, output):
        return management.Dissolve(source, output)

    def PairwiseClip(self, source, overlay, output):
        clip = _projected(overlay)
        return _save(_run('native:clip', INPUT=_projected(source, clip.crs()), OVERLAY=clip), output)

    def PairwiseIntersect(self, sources, output, join_attributes='ALL'):
        if len(sources) != 2 or join_attributes != 'ALL':
            raise NotImplementedError('LETO intersection requires two layers and all attributes.')
        left = _projected(sources[0])
        right = _projected(sources[1], left.crs())
        # Provider FIDs are not scientific attributes; avoid duplicate fid fields.
        right_fields = [f.name() for f in right.fields() if f.name().lower() != 'fid']
        return _save(_run('native:intersection', INPUT=left, OVERLAY=right,
                          INPUT_FIELDS=[f.name() for f in left.fields() if f.name().lower() != 'fid'],
                          OVERLAY_FIELDS=right_fields, OVERLAY_FIELDS_PREFIX='overlay_'), output)

    def PairwiseErase(self, source, overlay, output):
        left = _projected(source)
        return _save(_run('native:difference', INPUT=left, OVERLAY=_projected(overlay, left.crs())), output)

    def PairwiseBuffer(self, source, output, distance, dissolve_option='NONE', method='PLANAR'):
        if method != 'PLANAR':
            raise NotImplementedError('Only projected PLANAR buffers have a validated implementation.')
        layer = _projected(source)
        if not math.isclose(SpatialReference(layer.crs()).metersPerUnit, 1.0):
            raise ValueError('Buffering requires the model CRS in meters.')
        field_names = {f.name() for f in layer.fields()}
        buffer_field = '_leto_buffer_m'
        memory = layer.materialize(q.QgsFeatureRequest())
        if not memory.startEditing():
            raise RuntimeError('Cannot create buffer distance field.')
        memory.addAttribute(q.QgsField(buffer_field, QVariant.Double))
        index = memory.fields().indexFromName(buffer_field)
        for feature in memory.getFeatures():
            text = feature[distance] if distance in field_names else distance
            memory.changeAttributeValue(feature.id(), index, _distance_meters(text))
        if not memory.commitChanges():
            raise RuntimeError('Cannot commit buffer distance field.')
        result = _run('native:buffer', INPUT=memory, DISTANCE=q.QgsProperty.fromField(buffer_field),
                      SEGMENTS=BUFFER_SEGMENTS, END_CAP_STYLE=0, JOIN_STYLE=0, MITER_LIMIT=2,
                      DISSOLVE=dissolve_option == 'ALL')
        return _save(result, output)

    def GenerateNearTable(self, in_features, near_features, out_table, search_radius='',
                          location='NO_LOCATION', angle='NO_ANGLE', closest='CLOSEST',
                          closest_count=0, method='PLANAR'):
        if method.upper() != 'PLANAR':
            raise NotImplementedError('GEODESIC nearest-boundary search has no proven equivalent '
                                      'in this port. Do not substitute planar distance. '
                                      'Use the ArcPy version for proximity_merge or GEODESIC imputation.')
        left, right = _vector(in_features), _vector(near_features)
        if left.crs() != right.crs() or left.crs().isGeographic():
            raise ValueError('PLANAR near search requires matching projected coordinate systems.')
        factor = SpatialReference(left.crs()).metersPerUnit
        radius = _distance_meters(search_radius, factor) / factor if search_radius else None
        targets = list(right.getSelectedFeatures() if str(near_features) in _aliases else right.getFeatures())
        geometries = {f.id(): f.geometry() for f in targets}
        index = q.QgsSpatialIndex(iter(targets), flags=q.QgsSpatialIndex.FlagStoreFeatureGeometries)
        same = left.source() == right.source()
        sources = left.getSelectedFeatures() if str(in_features) in _aliases else left.getFeatures()
        rows, ties = [], []
        for feature in sources:
            geom = feature.geometry()
            if radius is not None:
                rect = geom.boundingBox()
                rect.grow(radius)
                candidates = index.intersects(rect)
            else:
                candidates = index.nearestNeighbor(geom, 2 if same else 1)
                candidates = [fid for fid in candidates if not (same and fid == feature.id())]
                if candidates:
                    distance = min(geom.distance(geometries[fid]) for fid in candidates)
                    rect = geom.boundingBox()
                    rect.grow(distance + max(1e-9, distance * 1e-12))
                    candidates = index.intersects(rect)
            pairs = sorted((geom.distance(geometries[fid]), fid) for fid in candidates
                           if not (same and fid == feature.id()))
            pairs = [(d, fid) for d, fid in pairs if d >= 0 and (radius is None or d <= radius)]
            if closest == 'CLOSEST' and pairs:
                tied = [(d, fid) for d, fid in pairs if math.isclose(d, pairs[0][0], rel_tol=1e-12, abs_tol=1e-9)]
                if len(tied) > 1:
                    ties.extend({'IN_FID': feature.id(), 'NEAR_FID': fid, 'NEAR_DIST': d} for d, fid in tied)
                pairs = pairs[:closest_count or 1]
            elif closest != 'ALL':
                if closest != 'CLOSEST':
                    raise ValueError(closest)
            rows.extend({'IN_FID': feature.id(), 'NEAR_FID': fid, 'NEAR_DIST': d} for d, fid in pairs)
        dest = _file(out_table, 'table')
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        if ties:
            tie_path = str(Path(dest).with_suffix('.ties.csv'))
            pd.DataFrame(ties).to_csv(tie_path, index=False)
            raise RuntimeError('Nearest-donor ties require baseline comparison; candidates written to '
                               f'{tie_path}. No arbitrary donor was selected.')
        pd.DataFrame(rows, columns=['IN_FID', 'NEAR_FID', 'NEAR_DIST']).to_csv(dest, index=False)
        return dest


def _distance_meters(value, map_unit_factor=1.0):
    if isinstance(value, (float, int)):
        return float(value) * map_unit_factor
    match = re.fullmatch(r'\s*([+-]?[\d.]+)\s*([A-Za-z]*)\s*', str(value))
    if not match:
        raise ValueError(f'Invalid distance: {value}')
    number, unit = float(match[1]), match[2].lower()
    factors = {'': map_unit_factor, 'meter': 1, 'meters': 1, 'm': 1,
               'foot': 0.3048, 'feet': 0.3048, 'ft': 0.3048,
               'kilometers': 1000, 'kilometer': 1000, 'km': 1000}
    if unit not in factors:
        raise ValueError(f'Unsupported distance unit: {unit}')
    return number * factors[unit]


management = _Management()
analysis = _Analysis()
conversion = _Conversion()
sa = _SpatialAnalyst()


def preflight(config, run_folder):
    """Fail before raster processing for unsupported spatial semantics."""
    import json
    initialize()
    if config.small_stands['remaining_small_stand_action'] == 'proximity_merge':
        raise NotImplementedError('proximity_merge requires GEODESIC near search. '
                                  'Keep this configuration on the original ArcPy backend until validated.')
    if str(config.inventory['nearest_method']).upper() != 'PLANAR':
        raise NotImplementedError('Only PLANAR inventory near search is supported by the QGIS port.')
    sr = Describe(config.treemap_raster).spatialReference
    if sr.type != 'Projected' or not math.isclose(sr.metersPerUnit, 1.0):
        raise ValueError('The baseline segmentation grid must be projected in meters.')
    # Verify original TreeMap donor attributes before doing expensive work.
    names = {field.name.upper() for field in ListFields(config.treemap_lookup_raster)}
    if 'VALUE' not in names or not names.intersection({'PLT_CN', 'PLTCN', 'PLOT_CN'}):
        raise ValueError('TreeMap lookup requires VALUE and PLT_CN raster attributes.')
    metadata = {'backend': 'PyQGIS/GDAL', 'qgis_version': q.Qgis.QGIS_VERSION,
                'gdal_version': gdal.VersionInfo('--version'), 'buffer_segments_per_quadrant': BUFFER_SEGMENTS,
                'source_grid': sr.exportToString(), 'polygonization_connectivity': 4,
                'categorical_resampling': 'nearest', 'polygon_mask': 'cell_center',
                'validation_status': 'Requires Maine ArcPy baseline comparison',
                'unsupported': ['GEODESIC near search', 'proximity_merge'],
                'nearest_tie_policy': 'fail and export candidates'}
    Path(run_folder, 'Backend_Settings.json').write_text(json.dumps(metadata, indent=2))
