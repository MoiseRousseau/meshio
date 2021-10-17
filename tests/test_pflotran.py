import pathlib

import numpy as np
import pytest

#import meshio
import sys
import os
path = os.getcwd() + '/../src/'
sys.path.insert(0,path)
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

@pytest.mark.parametrize(
    "filename, cell_set_name, cell_sets, n_bcs",
    [
        ("mixed_region.h5",
             "tetra", #name of the cell sets
             {"tetra":3}, 
             {"triangles":12, "quad":15},
         ),
    ]
)
def test_io_region_h5(filename, cell_set_name, cell_sets, n_bcs):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "pflotran" / filename
    mesh = meshio.read(filename)
    print(mesh.cell_sets[cell_set_name])
    cell_set_to_test = mesh.cell_sets[cell_set_name]
    for elem_type, n_elem in cell_set_to_test:
        if elem_type in cell_sets.keys():
            assert len(n_elem) == cell_sets[elem_type]
        else:
            assert len(n_elem) == 0
    for elem_type, n_elem in n_bcs.items():
        for x in mesh.cells:
           if x[0] == elem_type:
               assert len(x[1]) == n_elem
    #TODO: output
    return
    
@pytest.mark.parametrize(
    "meshname, regnames, refmesh",
    [
        ("mixed.ugi",
             ["tetra.vs", "bc.ss"], #name of the cell sets
             "mixed_region.h5", 
         ),
    ]
)
def test_io_region_ascii(meshname, regnames, refmesh):
    #TODO: output
    this_dir = pathlib.Path(__file__).resolve().parent
    refmesh = this_dir / "meshes" / "pflotran" / refmesh
    reference = meshio.read(refmesh)
    meshname = this_dir / "meshes" / "pflotran" / meshname
    mesh = meshio.read(meshname)
    for region in regnames:
        region = this_dir / "meshes" / "pflotran" / region
        meshio.pflotran.read_ascii_group(region, mesh)
    for elem_type, elem in reference.cells:
        for elem_type2, elem2 in mesh.cells:
            if elem_type == elem_type2:
                assert np.all(np.isclose(elem, elem2))
                break
    for cell_set_name, cell_set in reference.cell_sets.items():
        cell_set2 = mesh.cell_sets[cell_set_name]
        for elem_type, elem in cell_set:
	          for elem_type2, elem2 in cell_set2:
	              if elem_type == elem_type2:
	                  assert np.all(np.isclose(elem, elem2))
	                  break
    return

def test_write_material_ids_h5():
    #TODO
    return

