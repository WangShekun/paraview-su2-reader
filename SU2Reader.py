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
# - Parsing is vectorized with numpy.fromstring for large meshes; a per-line
#   fallback preserves precise error messages for malformed input.
# -----------------------------------------------------------------------------

from paraview.util.vtkAlgorithm import (
    smproxy, smproperty, smdomain, smhint, VTKPythonAlgorithmBase
)

import vtk
from vtk.util import numpy_support
import numpy as np


# vtkIdType size (4 or 8 bytes) depends on how VTK was built.
# Match the numpy dtype to avoid an extra copy/conversion in numpy_to_vtkIdTypeArray.
_ID_DTYPE = np.int64 if vtk.vtkIdTypeArray().GetDataTypeSize() == 8 else np.int32


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
# -----------------------------------------------------------------------------
def _strip_arrays(ds):
    if ds is None:
        return
    ds.GetPointData().Initialize()
    ds.GetCellData().Initialize()
    ds.GetFieldData().Initialize()


# -----------------------------------------------------------------------------
# Line / KV helpers
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


# -----------------------------------------------------------------------------
# Per-line fallback parsers (used when bulk numpy parse cannot reshape uniformly)
# -----------------------------------------------------------------------------
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
    if etype not in SU2_NNODES:
        raise RuntimeError(f"Unsupported {context} type {etype}")
    if allowed_types is not None and etype not in allowed_types:
        raise RuntimeError(f"Unsupported {context} type {etype}")

    expected = 1 + SU2_NNODES[etype]
    if len(parts) < expected:
        raise RuntimeError(
            f"{context} type {etype} has {len(parts) - 1} nodes, expected {SU2_NNODES[etype]}"
        )

    return etype, parts[1:expected]


# -----------------------------------------------------------------------------
# Bulk parsers (numpy-vectorized for big meshes)
# -----------------------------------------------------------------------------
def _bulk_parse_points(lines, start, npoin, ndime):
    """Parse npoin lines starting at lines[start] into a (npoin, 3) float64 array."""
    if npoin == 0:
        return np.zeros((0, 3), dtype=np.float64)
    if start + npoin > len(lines):
        raise RuntimeError(
            f"Expected {npoin} point lines, only {len(lines) - start} available"
        )

    pts = np.zeros((npoin, 3), dtype=np.float64)
    block = lines[start:start + npoin]

    # Fast path: uniform column count -> single np.fromstring on a joined buffer
    first_ncols = len(block[0].split())
    if first_ncols >= ndime:
        try:
            arr = np.fromstring(" ".join(block), sep=" ", dtype=np.float64)
        except Exception:
            arr = None
        if arr is not None and arr.size == first_ncols * npoin:
            pts[:, :ndime] = arr.reshape(npoin, first_ncols)[:, :ndime]
            return pts

    # Fallback: per-line (preserves precise error messages on malformed rows)
    for k, line in enumerate(block):
        pts[k, :ndime] = _parse_point_line(line, ndime, k)
    return pts


def _bulk_parse_elements(lines, start, nelem, allowed_types, context):
    """Parse nelem element lines into a list of (su2_etype, conns_np) groups.

    `conns_np` is a contiguous (count, nnodes) int64 array. Groups are returned
    in first-encounter order. For a homogeneous block (one element type, same
    column count on every row) this returns a single group via a single
    np.fromstring call -- the fast path used by typical CFD meshes.
    """
    if nelem == 0:
        return []
    if start + nelem > len(lines):
        raise RuntimeError(
            f"Expected {nelem} {context} lines, only {len(lines) - start} available"
        )

    block = lines[start:start + nelem]

    # Fast path: assume homogeneous (uniform column count + single etype)
    first_ncols = len(block[0].split())
    if first_ncols >= 2:
        try:
            arr = np.fromstring(" ".join(block), sep=" ", dtype=np.int64)
        except Exception:
            arr = None
        if arr is not None and arr.size == first_ncols * nelem:
            arr2 = arr.reshape(nelem, first_ncols)
            etypes = arr2[:, 0]
            if int(etypes.min()) == int(etypes.max()):
                etype = int(etypes[0])
                if etype not in SU2_NNODES:
                    raise RuntimeError(f"Unsupported {context} type {etype}")
                if allowed_types is not None and etype not in allowed_types:
                    raise RuntimeError(f"Unsupported {context} type {etype}")
                n = SU2_NNODES[etype]
                if first_ncols < 1 + n:
                    raise RuntimeError(
                        f"{context} type {etype} has {first_ncols - 1} nodes, expected {n}"
                    )
                conns = np.ascontiguousarray(arr2[:, 1:1 + n])
                return [(etype, conns)]

    # Fallback: per-line, group by etype preserving first-encounter order
    type_buf = {}
    order = []
    for k, line in enumerate(block):
        etype, conn = _parse_element_line(line, f"{context} {k}", allowed_types)
        if etype not in type_buf:
            type_buf[etype] = []
            order.append(etype)
        type_buf[etype].append(conn)
    return [(et, np.asarray(type_buf[et], dtype=np.int64)) for et in order]


# -----------------------------------------------------------------------------
# Core SU2 reader
# -----------------------------------------------------------------------------
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
    # Tolerate extended forms like 'NPOIN= 1234 5678' (only the count matters)
    npoin = int(_parse_kv_value(lines[i], "NPOIN").split()[0])
    pts = _bulk_parse_points(lines, i + 1, npoin, ndime)

    # ---- NELEM ----
    i_elem = find("NELEM")
    if i_elem is None:
        raise RuntimeError("NELEM not found")
    nelem = int(_parse_kv_value(lines[i_elem], "NELEM"))
    vol_groups = _bulk_parse_elements(
        lines, i_elem + 1, nelem, set(VTK_TYPE_MAP.keys()), "volume element"
    )

    # ---- NMARK ----
    i = find("NMARK")
    if i is None:
        raise RuntimeError("NMARK not found")
    nmark = int(_parse_kv_value(lines[i], "NMARK"))

    marker_groups = []
    p = i + 1
    allowed_marker_types = {3} if ndime == 2 else {5, 9}
    for _ in range(nmark):
        tag = _parse_kv_value(_require_data_line(lines, p, "MARKER_TAG"), "MARKER_TAG")
        p += 1
        ne = _parse_count_line(lines, p, "MARKER_ELEMS")
        p += 1
        groups = _bulk_parse_elements(
            lines, p, ne, allowed_marker_types, f"marker '{tag}' element"
        )
        p += ne
        marker_groups.append((tag, groups))

    return ndime, pts, vol_groups, marker_groups


# -----------------------------------------------------------------------------
# Build vtkUnstructuredGrid from grouped (numpy) connectivity
# -----------------------------------------------------------------------------
def build_ugrid(points, groups):
    """points: (N, 3) float64; groups: list of (su2_etype, (n_e, k) int array)."""
    ug = vtk.vtkUnstructuredGrid()

    # ---- Points: bulk insert ----
    pts_arr = np.ascontiguousarray(points, dtype=np.float64)
    vtk_pts = vtk.vtkPoints()
    vtk_pts.SetData(numpy_support.numpy_to_vtk(pts_arr, deep=True))
    ug.SetPoints(vtk_pts)

    # ---- Pre-pass: vectorized validation + size totals ----
    npts = points.shape[0]
    n_cells = 0
    total_ids = 0
    for etype, conns in groups:
        n_e, k = conns.shape
        if n_e == 0:
            continue
        cmin = int(conns.min())
        cmax = int(conns.max())
        if cmin < 0 or cmax >= npts:
            raise RuntimeError(
                f"Node index out of range for cell type {etype}; "
                f"expected 0 <= nid < {npts}, got [{cmin}, {cmax}]"
            )
        n_cells += n_e
        total_ids += n_e * (k + 1)

    if n_cells == 0:
        _strip_arrays(ug)
        return ug

    # ---- Build flat cell connectivity [k, id0, ..., id(k-1), k, id0, ...] ----
    flat = np.empty(total_ids, dtype=_ID_DTYPE)
    ctypes_arr = np.empty(n_cells, dtype=np.uint8)

    flat_off = 0
    cell_off = 0
    for etype, conns in groups:
        n_e, k = conns.shape
        if n_e == 0:
            continue
        size = n_e * (k + 1)
        view = flat[flat_off:flat_off + size].reshape(n_e, k + 1)
        view[:, 0] = k
        view[:, 1:] = conns
        flat_off += size

        ctypes_arr[cell_off:cell_off + n_e] = VTK_TYPE_MAP[etype]
        cell_off += n_e

    id_arr = numpy_support.numpy_to_vtkIdTypeArray(flat, deep=True)
    cells = vtk.vtkCellArray()
    cells.SetCells(n_cells, id_arr)

    vtk_ctypes = numpy_support.numpy_to_vtk(
        ctypes_arr, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR
    )
    ug.SetCells(vtk_ctypes, cells)

    _strip_arrays(ug)
    return ug


# -----------------------------------------------------------------------------
# Build a compacted boundary/surface grid for one marker
# (point set is reduced to nodes actually referenced by the marker)
# -----------------------------------------------------------------------------
def build_compacted(global_pts, groups):
    if not groups:
        return None
    groups = [(et, c) for et, c in groups if c.size]
    if not groups:
        return None

    # Union of referenced node IDs (sorted) -> remap via searchsorted (O(N log N))
    all_ids = np.unique(np.concatenate([c.ravel() for _, c in groups]))
    if all_ids.size:
        npts = global_pts.shape[0]
        if int(all_ids.min()) < 0 or int(all_ids.max()) >= npts:
            raise RuntimeError(
                f"Marker node index out of range; expected 0 <= nid < {npts}"
            )

    local_pts = global_pts[all_ids, :]
    new_groups = [(et, np.searchsorted(all_ids, c)) for et, c in groups]
    return build_ugrid(local_pts, new_groups)


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

        ndime, pts, vol_groups, marker_groups = read_su2(self._filename)

        # ---- Output 0: Volume mesh ----
        out0 = vtk.vtkUnstructuredGrid.GetData(outInfo, 0)
        out0.ShallowCopy(build_ugrid(pts, vol_groups))

        # ---- Output 1: Boundary meshes (MultiBlock) ----
        mb = vtk.vtkMultiBlockDataSet()
        blocks = []
        for tag, groups in marker_groups:
            surf = build_compacted(pts, groups)
            if surf is not None:
                blocks.append((tag, surf))

        mb.SetNumberOfBlocks(len(blocks))
        for i, (tag, b) in enumerate(blocks):
            mb.SetBlock(i, b)
            mb.GetMetaData(i).Set(vtk.vtkCompositeDataSet.NAME(), tag)

        out1 = vtk.vtkMultiBlockDataSet.GetData(outInfo, 1)
        out1.ShallowCopy(mb)

        return 1
