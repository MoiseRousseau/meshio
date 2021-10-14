import pathlib

import numpy as np
import pytest

import meshio

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
        helpers.hex_mesh,
        helpers.polyhedron_mesh,
    ],
)
def test_pflotran(mesh):
    def writer(*args, **kwargs):
        return meshio.pflotran.write(*args, **kwargs)

    helpers.write_read(writer, meshio.pflotran.read, mesh, 1.0e-12)


@pytest.mark.parametrize(
    "filename, sum_coordinate, ref_num_cells",
    [
        ("mixed.ugi", [48.75,66.25,61.25], {
            "tetra":3, 
            "pyramid":6, 
            "wedge":3,
            "hexahedron":3}),
        ("mixed.h5", [48.75,66.25,61.25], {
            "tetra":3, 
            "pyramid":6, 
            "wedge":3,
            "hexahedron":3}),
    ],
)
def test_reference_file(filename, sum_coordinate, ref_num_cells):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "pflotran" / filename

    mesh = meshio.read(filename)
    s = np.sum(mesh.points,axis=0)
    assert abs(s[0] - ref_sum[0]) < 1e-6
    assert abs(s[1] - ref_sum[1]) < 1e-6
    assert abs(s[2] - ref_sum[2]) < 1e-6
    for elem_type, num in ref_num_cells.items:
        assert len(mesh.get_cells_type(elem_type)) == num

def test_region():
    #TODO
    return

def test_material_ids():
    #TODO
    return

def test_bc():
    #TODO
    return
