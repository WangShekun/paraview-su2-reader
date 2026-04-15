# SU2Reader.py
# -----------------------------------------------------------------------------
# ParaView Python Reader Plugin for SU2 mesh
#
# Outputs:
#   Output 0: Volume Mesh        (vtkUnstructuredGrid)
#   Output 1: Boundary Meshes    (vtkMultiBlockDataSet, one block per marker)
#
# Design notes:
# - Volume and marker connectivity follow SU2's VTK-style element sizes
# - Trailing legacy IDs are ignored when present
# - Marker names are kept for MultiBlock inspection
# - No data arrays are attached (pure geometry)
# -----------------------------------------------------------------------------

from paraview.util.vtkAlgorithm import (
    smproxy, smproperty, smdomain, smhint, VTKPythonAlgorithmBase
)

import vtk
import numpy as np


# -----------------------------------------------------------------------------
# SU2 element type -> VTK cell type
# -----------------------------------------------------------------------------
VTK_TYPE_MAP = {
    3: 3,    # line
    5: 5,    # triangle
    9: 9,    # quad
    10: 10,  # tet
    12: 12,  # hex
    13: 13,  # prism (wedge)
    14: 14   # pyramid
}

# -----------------------------------------------------------------------------
# Expected node counts for SU2 elements
# (used to read only the connectivity nodes and ignore trailing legacy IDs)
# -----------------------------------------------------------------------------
SU2_NNODES = {
    3: 2,
    5: 3,
    9: 4,
    10: 4,
    12: 8,
    13: 6,
    14: 5
}


# -----------------------------------------------------------------------------
# Utility: remove all data arrays (geometry only)
# NOTE:
# - This prevents scalar coloring from data arrays
# - It does NOT control ParaView's MultiBlock representation behavior
# -----------------------------------------------------------------------------
def _strip_arrays(ds):
    if ds is None:
        return
    ds.GetPointData().Initialize()
    ds.GetCellData().Initialize()
    ds.GetFieldData().Initialize()


# -----------------------------------------------------------------------------
# Core SU2 reader
# -----------------------------------------------------------------------------
def _clean_line(line):
    return line.split("%", 1)[0].strip()


def _parse_kv_value(line, key):
    head, sep, tail = line.partition("=")
    if sep != "=" or head.strip() != key:
        raise RuntimeError(f"Expected '{key}=' line, got: {line}")
    value = tail.strip()
    if not value:
        raise RuntimeError(f"Missing value for {key}")
    return value


def _parse_count_line(lines, index, key):
    if index >= len(lines):
        raise RuntimeError(f"Expected '{key}=' line, got EOF")
    return int(_parse_kv_value(lines[index], key))


def _require_data_line(lines, index, context):
    if index >= len(lines):
        raise RuntimeError(f"Unexpected EOF while reading {context}")
    return lines[index]


def _parse_point_line(line, ndime, point_index):
    parts = line.split()
    if len(parts) < ndime:
        raise RuntimeError(
            f"Point {point_index} has {len(parts)} values, expected at least {ndime}"
        )
    return [float(value) for value in parts[:ndime]]


def _parse_element_line(line, context, allowed_types=None):
    parts = [int(value) for value in line.split()]
    if not parts:
        raise RuntimeError(f"Empty {context} line")

    etype = parts[0]
    if allowed_types is not None and etype not in allowed_types:
        raise RuntimeError(f"Unsupported {context} type {etype}")
    if etype not in SU2_NNODES:
        raise RuntimeError(f"Unsupported {context} type {etype}")

    expected = 1 + SU2_NNODES[etype]
    if len(parts) < expected:
        raise RuntimeError(
            f"{context} type {etype} has {len(parts) - 1} nodes, expected {SU2_NNODES[etype]}"
        )

    return etype, parts[1:expected]


def _validate_conn(conn, npoints, context):
    if not conn:
        raise RuntimeError(f"Empty connectivity for {context}")
    for nid in conn:
        if nid < 0 or nid >= npoints:
            raise RuntimeError(
                f"Node index {nid} out of range for {context}; expected 0 <= nid < {npoints}"
            )


def read_su2(path):
    with open(path, "r") as f:
        lines = []
        for raw in f:
            line = _clean_line(raw)
            if line:
                lines.append(line)

    def find(key):
        for i, l in enumerate(lines):
            if l.startswith(key):
                return i
        return None

    # ---- NDIME ----
    i = find("NDIME")
    if i is None:
        raise RuntimeError("NDIME not found")
    ndime = int(_parse_kv_value(lines[i], "NDIME"))
    if ndime not in (2, 3):
        raise RuntimeError(f"Unsupported NDIME={ndime}; expected 2 or 3")

    # ---- NPOIN ----
    i = find("NPOIN")
    if i is None:
        raise RuntimeError("NPOIN not found")
    npoin = int(_parse_kv_value(lines[i], "NPOIN"))

    pts = np.zeros((npoin, 3))
    for k in range(npoin):
        line = _require_data_line(lines, i + 1 + k, f"point {k}")
        xyz = _parse_point_line(line, ndime, k)
        pts[k, :ndime] = xyz

    # ---- NELEM ----
    vol_conns, vol_ctypes = [], []
    i_elem = find("NELEM")
    if i_elem is None:
        raise RuntimeError("NELEM not found")
    nelem = int(_parse_kv_value(lines[i_elem], "NELEM"))
    for k in range(nelem):
        line = _require_data_line(lines, i_elem + 1 + k, f"volume element {k}")
        etype, conn = _parse_element_line(line, f"volume element {k}", VTK_TYPE_MAP)
        vol_conns.append(conn)
        vol_ctypes.append(VTK_TYPE_MAP[etype])

    # ---- NMARK ----
    i = find("NMARK")
    if i is None:
        raise RuntimeError("NMARK not found")
    nmark = int(_parse_kv_value(lines[i], "NMARK"))

    markers = {}
    p = i + 1
    allowed_marker_types = {3} if ndime == 2 else {5, 9}
    for _ in range(nmark):
        tag = _parse_kv_value(_require_data_line(lines, p, "MARKER_TAG"), "MARKER_TAG")
        p += 1
        ne = _parse_count_line(lines, p, "MARKER_ELEMS")
        p += 1

        elems = []
        for elem_idx in range(ne):
            line = _require_data_line(lines, p, f"marker element {elem_idx} for '{tag}'")
            etype, conn = _parse_element_line(
                line,
                f"marker element {elem_idx} for '{tag}'",
                allowed_marker_types,
            )
            elems.append([etype] + conn)
            p += 1

        markers[tag] = elems

    return ndime, pts, vol_conns, vol_ctypes, markers


# -----------------------------------------------------------------------------
# Build vtkUnstructuredGrid
# -----------------------------------------------------------------------------
def build_ugrid(points, conns, ctypes):
    ug = vtk.vtkUnstructuredGrid()

    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetNumberOfPoints(points.shape[0])
    for i, p in enumerate(points):
        vtk_pts.SetPoint(i, float(p[0]), float(p[1]), float(p[2]))
    ug.SetPoints(vtk_pts)

    for conn, ct in zip(conns, ctypes):
        _validate_conn(conn, points.shape[0], f"cell type {ct}")
        ids = vtk.vtkIdList()
        ids.SetNumberOfIds(len(conn))
        for j, nid in enumerate(conn):
            ids.SetId(j, int(nid))
        ug.InsertNextCell(int(ct), ids)

    _strip_arrays(ug)
    return ug


# -----------------------------------------------------------------------------
# Build compacted surface grid for one marker (3D markers: tri/quad)
# -----------------------------------------------------------------------------
def build_surface_compacted(global_pts, elems):
    used = set()
    conns, ctypes = [], []

    for e in elems:
        if not e:
            continue
        etype = int(e[0])
        if etype in (5, 9):  # tri / quad only
            conn = e[1:]
            _validate_conn(conn, global_pts.shape[0], f"marker surface type {etype}")
            conns.append(conn)
            ctypes.append(VTK_TYPE_MAP[etype])
            used.update(conn)

    if not conns:
        return None

    used = sorted(used)
    old2new = {o: i for i, o in enumerate(used)}

    local_pts = global_pts[used, :]
    local_conns = [[old2new[n] for n in c] for c in conns]

    return build_ugrid(local_pts, local_conns, ctypes)


# -----------------------------------------------------------------------------
# Build compacted boundary grid for one marker (2D only: LINE elements)
# -----------------------------------------------------------------------------
def build_boundary_compacted_2d(global_pts, elems):
    used = set()
    conns, ctypes = [], []

    for e in elems:
        if not e:
            continue
        if int(e[0]) == 3:  # line
            conn = e[1:3]
            if len(conn) == 2:
                _validate_conn(conn, global_pts.shape[0], "2D boundary line")
                conns.append(conn)
                ctypes.append(VTK_TYPE_MAP[3])
                used.update(conn)

    if not conns:
        return None

    used = sorted(used)
    old2new = {o: i for i, o in enumerate(used)}

    local_pts = global_pts[used, :]
    local_conns = [[old2new[n] for n in c] for c in conns]

    return build_ugrid(local_pts, local_conns, ctypes)


# -----------------------------------------------------------------------------
# ParaView Reader
# -----------------------------------------------------------------------------
@smproxy.reader(
    name="SU2 Mesh Reader",
    extensions="su2",
    file_description="SU2 Mesh"
)
@smhint.xml("""
  <Hints>
    <OutputPort index="0" name="Volume Mesh"/>
    <OutputPort index="1" name="Boundary Meshes"/>
    <NoDataArrays/>
  </Hints>
""")
class SU2MeshReader(VTKPythonAlgorithmBase):

    def __init__(self):
        super().__init__(nInputPorts=0, nOutputPorts=2)
        self._filename = None

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="su2", file_description="SU2 mesh")
    def SetFileName(self, fname):
        if fname != self._filename:
            self._filename = fname
            self.Modified()

    def FillOutputPortInformation(self, port, info):
        if port == 0:
            info.Set(vtk.vtkDataObject.DATA_TYPE_NAME(), "vtkUnstructuredGrid")
        else:
            info.Set(vtk.vtkDataObject.DATA_TYPE_NAME(), "vtkMultiBlockDataSet")
        return 1

    def RequestData(self, request, inInfo, outInfo):
        if not self._filename:
            return 0

        ndime, pts, conns, ctypes, markers = read_su2(self._filename)

        # ---- Output 0: Volume mesh ----
        out0 = vtk.vtkUnstructuredGrid.GetData(outInfo, 0)
        out0.ShallowCopy(build_ugrid(pts, conns, ctypes))

        # ---- Output 1: Boundary meshes (MultiBlock) ----
        mb = vtk.vtkMultiBlockDataSet()
        blocks = []

        for tag, elems in markers.items():
            if ndime == 2:
                surf = build_boundary_compacted_2d(pts, elems)
            else:
                surf = build_surface_compacted(pts, elems)

            if surf is not None:
                blocks.append((tag, surf))

        mb.SetNumberOfBlocks(len(blocks))
        for i, (tag, b) in enumerate(blocks):
            mb.SetBlock(i, b)
            # Keep marker names for ParaView inspection
            mb.GetMetaData(i).Set(
                vtk.vtkCompositeDataSet.NAME(), tag
            )

        out1 = vtk.vtkMultiBlockDataSet.GetData(outInfo, 1)
        out1.ShallowCopy(mb)

        return 1
