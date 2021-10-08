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
  for grp in s_grps:
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
  if "Regions" in mesh.keys(): 
    for reg in list(mesh["Regions"].keys()):
      if "Vertex Ids" in mesh[f"Regions/{reg}"].keys(): #face group
        grps[reg] = np.array(mesh[f"Regions/{reg}/Vertex Ids"])
      else:
        grps[reg] = np.array(mesh[f"Regions/{reg}/Cell Ids"])
  mesh.close()
  return vert, elem, grps


def read(filename):
    pft_to_meshio_order
    #read mesh and groups
    if filename.split('.') == "h5":
        vertices, elements, groups = read_h5(filename)
    else:
        vertices, elements, groups = read_ugi(filename)
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
    #create cell set
    
    #build mesh
    m = Mesh(
        points=vertices,
        cells=cells,
        cell_set=cell_data,
    )
    return



def write_h5(filename, mesh):
    #write vertices
    out = h5py.File(filename, 'w')
    out.create_dataset("Domain/Vertices", data=mesh.points)
    #write elements
    n_cells = 0
    for elem_type,cells in mesh.cells:
      if elem_type not in ["tetra", "pyramid", "wedge", "hexahedron"]:
        print(f"Element type {elem_type} not supported by PFLOTRAN, stop")
        return
      n_cells += len(cells)
    elements = np.zeros((n_cells,9), dtype='i8')
    count = 0
    for elem_type,cells in mesh.cells:
      elements[count:len(cells),0] = len(cells[0])
      elements[count:len(cells),1:len(cells[0])] = cells
      count += len(cells)
    out.create_dataset("Domain/Cells", data=elements)
    #write group
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



def write(filename, mesh):
    if filename.split('.')[-1] == "h5":
        write_h5(filename, mesh)
    else:
        write_ugi(filename, mesh)
    return


register("pflotran", [".h5", ".ugi"], read, {"pflotran": write})
