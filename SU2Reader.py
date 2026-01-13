# SU2Reader.py
# -----------------------------------------------------------------------------
# ParaView Python Reader Plugin for SU2 mesh
#
# Outputs:
#   Output 0: Volume Mesh        (vtkUnstructuredGrid)
#   Output 1: Boundary Meshes    (vtkMultiBlockDataSet, one block per marker)
#
# Design notes:
# - 3D import logic is strictly preserved
# - 2D support is added silently (no behavior change for 3D)
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
# (ONLY used to fix 2D meshes with trailing element IDs)
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
def read_su2(path):
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("%")]

    def find(key):
        for i, l in enumerate(lines):
            if l.startswith(key):
                return i
        return None

    # ---- NDIME ----
    i = find("NDIME")
    if i is None:
        raise RuntimeError("NDIME not found")
    ndime = int(lines[i].split("=")[1])

    # ---- NPOIN ----
    i = find("NPOIN")
    if i is None:
        raise RuntimeError("NPOIN not found")
    npoin = int(lines[i].split("=")[1])

    pts = np.zeros((npoin, 3))
    for k in range(npoin):
        xyz = list(map(float, lines[i + 1 + k].split()[:ndime]))
        pts[k, :ndime] = xyz

    # ---- NELEM ----
    vol_conns, vol_ctypes = [], []
    i_elem = find("NELEM")
    if i_elem is not None:
        nelem = int(lines[i_elem].split("=")[1])
        for k in range(nelem):
            parts = list(map(int, lines[i_elem + 1 + k].split()))
            etype = parts[0]

            # Ignore unsupported element types defensively
            if etype not in VTK_TYPE_MAP:
                continue

            # 2D fix: trim possible trailing element IDs
            # 3D path: EXACTLY identical to original behavior
            if ndime == 2 and etype in SU2_NNODES:
                conn = parts[1:1 + SU2_NNODES[etype]]
            else:
                conn = parts[1:]

            vol_conns.append(conn)
            vol_ctypes.append(VTK_TYPE_MAP[etype])

    # ---- NMARK ----
    i = find("NMARK")
    if i is None:
        raise RuntimeError("NMARK not found")
    nmark = int(lines[i].split("=")[1])

    markers = {}
    p = i + 1
    for _ in range(nmark):
        tag = lines[p].split("=")[1].strip()
        p += 1
        ne = int(lines[p].split("=")[1])
        p += 1

        elems = []
        for _ in range(ne):
            e = list(map(int, lines[p].split()))

            # 2D fix: marker elements may also carry trailing IDs
            if ndime == 2 and e and e[0] in SU2_NNODES:
                e = e[:1 + SU2_NNODES[e[0]]]

            elems.append(e)
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
        ids = vtk.vtkIdList()
        ids.SetNumberOfIds(len(conn))
        for j, nid in enumerate(conn):
            ids.SetId(j, int(nid))
        ug.InsertNextCell(int(ct), ids)

    _strip_arrays(ug)
    return ug


# -----------------------------------------------------------------------------
# Build compacted surface grid for one marker (unchanged 3D behavior)
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
                vtk.vtkCompositeDataSet.NAME(), f"marker:{tag}"
            )

        out1 = vtk.vtkMultiBlockDataSet.GetData(outInfo, 1)
        out1.ShallowCopy(mb)

        return 1

