"""
I/O for PFLOTRAN HDF5 implicit unstructured grid format.
"""
import datetime
import logging

import numpy as np
import h5py

from ..__about__ import __version__
from .._exceptions import WriteError
from .._files import open_file
from .._helpers import register
from .._mesh import CellBlock, Mesh


def read_ugi(src):
  mesh = open(src,'r')
  n_e, n_v = [int(x) for x in mesh.readline().split()]
  elem = np.zeros((n_e,9),dtype='i8')
  for iline in range(n_e):
    line = [int(x) for x in mesh.readline().split()[1:]]
    elem[iline,0] = len(line)
    elem[iline,1:len(line)+1] = line
  vert = np.zeros((n_v,3),dtype='i8')
  for iline in range(n_v):
    vert[iline] = [float(x) for x in mesh.readline().split()]
  mesh.close()
  #groups
  folder = '/'.join(src.split('/')[:-1])
  list_file = os.listdir(folder)
  s_grps = [x for x in list_file if x.split('.')[-1] == "ss"]
  v_grps = [x for x in list_file if x.split('.')[-1] == "vs"]
  grps = {}
  data = {}
  for grp in s_grps: #TODO: create new elements
    name = grp.split("/")[-1].split(".")[:-1][0]
    src = open(grp, 'r')
    n_face = int(src.readline())
    array = np.zeros((n_face,4), dtype='i8')
    for i in range(n_face):
      line = src.readline().split()
      array[i] = [int(x) for x in line[1:]]
    grps[name] = array
    src.close()
  for grp in v_grps:
    name = grp.split("/")[-1].split(".")[:-1][0]
    grps[name] = np.genfromtxt(grp, skip_header=1)
  return vert, elem, grps
  
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
  if "Materials" in mesh.keys(): 
      data["Materials"] = np.array(mesh[f"Regions/Materials/Material Ids"])
      #TODO: this suppose this is defined on all the cells numbered for 0 to N-1...
  mesh.close()
  return vert, elem, grps, data


def read(filename):
    #read mesh and groups
    ext = filename.split('.')
    if ext == "h5":
        vertices, elements, groups, data = read_h5(filename)
    elif ext == "ugi":
        vertices, elements, groups, data = read_ugi(filename)
    else:
        print("Meshio cannot read PFLOTRAN explicit mesh")
        return
    #convert to meshio format
    tets = elements[elements[0] == 4][:,1:5]
    pyrs = elements[elements[0] == 5][:,1:6]
    wedges = elements[elements[0] == 6][:,1:7]
    hexs = elements[elements[0] == 8][:,1:]
    cells = [("tetra", tets),
             ("pyramid", pyrs),
             ("wedge", wedges),
             ("hexahedron", hexs),
    ]
    #correspondance for groups
    index = np.arange(len(elements))
    tets = index[elements[0] == 4][:,1:5]
    pyrs = index[elements[0] == 5][:,1:6]
    wedges = index[elements[0] == 6][:,1:7]
    hexs = index[elements[0] == 8][:,1:]
    indexes = [("tetra", tets),
             ("pyramid", pyrs),
             ("wedge", wedges),
             ("hexahedron", hexs),
    ]
    #cells set
    cell_sets = {}
    for grp_name, grp_cells in grps.items():
        if len(grp_cells.shape) == 2:
            continue #boundary connections
        grp_cell_set = []
        for (type_name, type_index) in indexes:
            grp_cell_set.append(
                (
                type_name, 
                np.argwhere(np.isin(type_index, grp_cells)) 
                )
            )
        cell_sets[grp_name] = grp_cell_set
    #boundary connections
    for bc_name, bc_conn in grps.items():
        if len(bc_conn.shape) == 1:
            continue #already done above
        tri = bc_conn[bc_conn[0] == 3][:,1:4]
        quad = bc_conn[bc_conn[0] == 4][:,1:5]
        cells.append( ("triangle",tri) )
        cells.append( ("quad",quad) )
    #build mesh
    m = Mesh(
        points=vertices,
        cells=cells,
        cell_sets=groups,
    )
    return



def write_h5(filename, mesh):
    #get elements
    n_cells = 0
    for elem_type,cells in mesh.cells:
        #TODO write with BC
        if elem_type not in ["triangle", "quad", "tetra", "pyramid", 
                             "wedge", "hexahedron"]:
            if "polyhedron" in elem_type:
                write_h5_explicit(filename, mesh)
                return
            print(f"Element type {elem_type} not supported by PFLOTRAN, stop")
            return
        n_cells += len(cells)
    elements = np.zeros((n_cells,9), dtype='i8')
    count = 0
    for elem_type,cells in mesh.cells:
        elements[count:len(cells),0] = len(cells[0])
        elements[count:len(cells),1:len(cells[0])] = cells
        count += len(cells)
    #write vertices and elements
    out = h5py.File(filename, 'w')
    out.create_dataset("Domain/Vertices", data=mesh.points)
    out.create_dataset("Domain/Cells", data=elements)
    #write group
    for grp_name, grp in mesh.cell_sets: #TODO review
        out.create_dataset(f"Regions/{grp_name}/Cell Ids", data=elements)
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
      code = elem_type[0].uppercase()
      for cell in cells
        out.write(code)
        for index in cell:
          out.write(f" {index}")
        out.write('\n')
    out.close()
    #write points
    for p in mesh.points:
      out.write(f"{p[0]} {p[1]} {p[2]}\n")
    return

def write_uge(filename, mesh):
    #TODO
    return

def write_h5_explicit(filename, mesh):
    #TODO
    return



def write(filename, mesh):
    ext = filename.split('.')[-1]
    if ext == "h5":
        write_h5(filename, mesh)
    elif ext == "ugi":
        write_ugi(filename, mesh)
    elif ext == "uge":
        write_uge(filename, mesh)
    return


register("pflotran", [".h5", ".ugi", ".uge"], read, {"pflotran": write})
