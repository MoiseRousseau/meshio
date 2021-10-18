"""
I/O for PFLOTRAN HDF5 implicit unstructured grid format.
"""
import datetime
import logging
import os

import numpy as np
import h5py

from ..__about__ import __version__
from .._exceptions import WriteError
from .._files import open_file
from .._helpers import register
from .._mesh import CellBlock, Mesh


def read_ugi(filename):
    mesh = open(filename,'r')
    n_e, n_v = [int(x) for x in mesh.readline().split()]
    elem = np.zeros((n_e,9),dtype='i8')
    for iline in range(n_e):
        line = [int(x) for x in mesh.readline().split()[1:]]
        elem[iline,0] = len(line)
        elem[iline,1:len(line)+1] = line
    vert = np.zeros((n_v,3),dtype='f8')
    for iline in range(n_v):
        vert[iline] = [float(x) for x in mesh.readline().split()]
    mesh.close()
    return vert, elem, None, None #no groups, no data


def read_ascii_group(filename, mesh):
    surf = False
    if str(filename).split('.')[-1] == "ss":
        surf = True
    name = str(filename).split("/")[-1].split(".")[:-1][0]
    if surf: 
        src = open(filename, 'r')
        n_face = int(src.readline())
        grp = np.zeros((n_face,5), dtype='i8')
        cell_type = {"T":3, "Q":4}
        for i in range(n_face):
          line = src.readline().split()
          grp[i,0] = cell_type[line[0]]
          grp[i,1:cell_type[line[0]]+1] = [int(x) for x in line[1:]]
        src.close()
        add_bc_to_mesh(name, grp, mesh)
    else:
        grp = np.genfromtxt(filename)
        add_group_to_mesh(name, grp, mesh)
    return grp

  
def read_h5(src):
    mesh = h5py.File(src, 'r')
    vert = np.array(mesh["Domain/Vertices"])
    elem = np.array(mesh["Domain/Cells"])
    grps = {}
    data = {}
    if "Regions" in mesh.keys(): 
        for reg in list(mesh["Regions"].keys()):
            if "Vertex Ids" in mesh[f"Regions/{reg}"].keys():
              grps[reg] = np.array(mesh[f"Regions/{reg}/Vertex Ids"])
            else:
              grps[reg] = np.array(mesh[f"Regions/{reg}/Cell Ids"])
        if "Materials" in list(mesh["Regions"].keys()): 
            data["Materials"] = np.array(mesh[f"Regions/Materials/Material Ids"])
            #TODO: this suppose this is defined on all the cells numbered for 0 to N-1...
    mesh.close()
    return vert, elem, grps, data


def read(filename):
    #read mesh and groups
    ext = str(filename).split('.')[-1]
    if ext == "h5":
        vertices, elements, groups, data = read_h5(filename)
    elif ext == "ugi":
        vertices, elements, groups, data = read_ugi(filename)
    else:
        print("Warning! Unknown PFLOTRAN mesh extension, try ugi")
        vertices, elements, groups, data = read_ugi(filename)
    #convert to meshio format
    elem_type_code = {"tetra":4, "pyramid":5, "wedge":6, "hexahedron":8}
    cells = []
    for elem_type, elem_code in elem_type_code.items():
       cond = (elements[:,0] == elem_code)
       cells_of_type = elements[cond][:,1:elem_code+1]
       if len(cells_of_type):
           cells.append((elem_type, cells_of_type-1)) #-1 pflotran one based
    #create a cell data to store natural id
    index = np.arange(len(elements))+1
    indexes = []
    for elem_type, elem_code in elem_type_code.items():
       cond = (elements[:,0] == elem_code)
       cells_of_type = index[cond]
       if np.sum(cond):
           indexes.append(cells_of_type)
    #build mesh
    m = Mesh(
        points=vertices,
        cells=cells,
        cell_data={"pflotran_indexes":indexes}
    )
    if groups:
        cell_sets = {}
        for grp_name, grp_cells in groups.items():
            if len(grp_cells.shape) == 2:
                add_bc_to_mesh(grp_name, grp_cells, m)
            else:
                add_group_to_mesh(grp_name, grp_cells, m)
    for x in m.cells:
        if x[0] == "triangle":
            m.cell_data["pflotran_indexes"].append([-1 for i in x[1]])
        elif x[0] == "quad":
            m.cell_data["pflotran_indexes"].append([-1 for i in x[1]])
    return m


def add_group_to_mesh(grp_name, grp_cells, m):
    grp_cell_set = []
    elem_type_code = {"tetra":4, "pyramid":5, "wedge":6, "hexahedron":8}
    for (type_name, type_index) in zip(elem_type_code.keys(),
                                       m.cell_data["pflotran_indexes"]):
        cond = np.isin(type_index, grp_cells)
        data = np.argwhere(cond)
        grp_cell_set.append((type_name, data))
    m.cell_sets[grp_name] = grp_cell_set
    return


def add_bc_to_mesh(bc_name, bc_conn, m):
    tri = bc_conn[bc_conn[:,0] == 3][:,1:4]
    quad = bc_conn[bc_conn[:,0] == 4][:,1:5]
    added = False
    if len(tri): 
        for elem_type, elem in m.cells:
            if elem_type == "triangle":
                elem += tri
                added = True
                break
        if not added:
            m.cells.append( ("triangle",tri) )
    added = False
    if len(quad):
        for elem_type, elem in m.cells:
            if elem_type == "quad":
                elem += quad
                added = True
                break
        if not added:
            m.cells.append( ("quad",tri) )
    return


def write_h5(filename, mesh):
    #get elements
    n_cells = 0
    for elem_type,cells in mesh.cells:
        if elem_type not in ["triangle", "quad", "tetra", "pyramid", 
                             "wedge", "hexahedron"]:
            print(f"Element type {elem_type} not supported by PFLOTRAN, stop")
            return
        if elem_type in ["triangle", "quad"]:
            continue
        n_cells += len(cells)
    elements = np.zeros((n_cells,9), dtype='i8')
    count = 0
    for elem_type,cells in mesh.cells:
        if elem_type in ["triangle", "quad"]:
            continue
        elements[count:len(cells),0] = len(cells[0])
        elements[count:len(cells),1:len(cells[0])+1] = cells
        count += len(cells)
    #write vertices and elements
    out = h5py.File(filename, 'w')
    out.create_dataset("Domain/Vertices", data=mesh.points)
    out.create_dataset("Domain/Cells", data=elements+1) #pflotran is one based
    #write group
    for grp_name, grp in mesh.cell_sets:
        if ("triangle" in grp.keys() or
            "quad" in grp.keys):
            continue
        count = 0
        grp_indexes = []
        for elem_type, elem in mesh.cells:
            if elem_type not in grp.keys():
                count += len(elem)
                continue
            grp_indexes += list(grp[elem_type]+1+count)
        out.create_dataset(f"Regions/{grp_name}/Cell Ids", data=grp_indexes)
    #write bc
    for name, cell_set in mesh.cell_sets:
        for elem_type, elem in cell_set:
            bc = np.zeros((len(cell),5), dtype='i8')
            bc[:,0] = len(elem[0])
            bc[:,1:] = elem+1
        out.write(f"Region/{name}/Vertex Ids", dataset=bc)
    out.close()
    return

def write_ugi(filename, mesh):
    #determine number of cells and points
    n_pts = len(mesh.points)
    n_cells = 0
    for elem_type,cells in mesh.cells:
        if elem_type not in ["tetra", "pyramid", "wedge", "hexahedron"]:
            print(f"Element type {elem_type} not supported by PFLOTRAN, stop")
            return
        n_cells += len(cells)
        #out cells
        out = open(filename, 'w')
        out.write(f"{n_cells} {n_pts}\n")
    #write cells
    for elem_type,cells in mesh.cells:
      code = elem_type[0].upper()
      for cell in cells:
          out.write(code)
          for index in cell:
              out.write(f" {index+1}")
          out.write('\n')
    #write points
    for p in mesh.points:
        out.write(f"{p[0]} {p[1]} {p[2]}\n")
    out.close()
    #write groups
    _filename = str(filename)
    name_ext = _filename.split("/")[-1]
    folder = _filename[:len(_filename)-len(name_ext)+1]
    for grp_name, grp in mesh.cell_sets:
        count = 0
        grp_indexes = []
        for elem_type, elem in mesh.cells:
            if elem_type not in grp.keys():
                count += len(elem)
                continue
            grp_indexes += list(grp[elem_type]+1+count)
        out = open(folder+name, 'w')
        for i in grp_indexes:
            out.write(f'{i}\n')
        out.close()
    #write bc
    for name, cell_set in mesh.cell_sets:
        for elem_type, elem in cell_set:
            bc = np.zeros((len(cell),5), dtype='i8')
            bc[:,0] = len(elem[0])
            bc[:,1:len(elem[0])+1] = elem+1
        out = open(folder+name, 'w')
        for face in bc:
           if face[0] == 3:
               out.write("T ")
               for i in range(3):
                   out.write(f"{face[i]} ")
               out.write("\n")
           elif face[0] == 4:
               out.write("Q ")
               for i in range(4):
                   out.write(f"{face[i]} ")
               out.write("\n")
        out.close()
    return



def write(filename, mesh):
    ext = str(filename).split('.')[-1]
    if ext == "h5":
        write_h5(filename, mesh)
    elif ext == "ugi":
        write_ugi(filename, mesh)
    else:
        print("Warning! Unknown PFLOTRAN mesh extension, export as ugi")
        write_ugi(filename, mesh)
    return


register("pflotran", [".h5", ".ugi"], read, {"pflotran": write})
