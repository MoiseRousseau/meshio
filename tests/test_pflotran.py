import pathlib

import numpy as np
import pytest

#import meshio
import sys
import os
path = os.getcwd() + '/../src/'
sys.path.insert(0,path)
print(sys.path)
import meshio

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
        helpers.hex_mesh,
        #helpers.polyhedron_mesh,
    ],
)
def test_io_ascii(mesh):
    #Test write and read function for consistency
    helpers.write_read(meshio.pflotran.write, 
                       meshio.pflotran.read, 
                       mesh, 
                       1.0e-12,
                       extension="ugi")


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
        helpers.hex_mesh,
        #helpers.polyhedron_mesh,
    ],
)
def test_io_h5(mesh):
    helpers.write_read(meshio.pflotran._pflotran.write, 
                       meshio.pflotran._pflotran.read, 
                       mesh, 
                       1.0e-12, 
                       extension="h5")


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
def test_read_reference_file(filename, sum_coordinate, ref_num_cells):
    #test HDF5 / ASCII read
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "pflotran" / filename

    mesh = meshio.read(filename)
    s = np.sum(mesh.points,axis=0)
    print(mesh.points)
    assert abs(s[0] - sum_coordinate[0]) < 1e-6
    assert abs(s[1] - sum_coordinate[1]) < 1e-6
    assert abs(s[2] - sum_coordinate[2]) < 1e-6
    for elem_type, num in ref_num_cells.items():
        assert len(mesh.get_cells_type(elem_type)) == num


def test_io_region_h5():
    #TODO
    return

def test_write_material_ids_h5():
    #TODO
    return

def test_bc_h5():
    #TODO
    return

def test_io_region_ascii():
    #TODO
    return

def test_bc_ascii():
    #TODO
    return
