from netgen.occ import *
from ngsolve import *

"""
James Elgy - 2022:
Dual bar example for Netgen OCC geometry mesh generation.
This example illustrates that lists can be used for multiple objects.
"""



# Setting mur, sigma, and defining the top level object name:
# For multiple objects then object_name, mur, and sigma should be given as parameter lists.
# These lists are used in the generation of the associated .geo files.
material_name = ['bar1', 'bar2']
mur = [1, 1]
sigma = [1e6, 1e8]
alpha = 0.001

# Setting Boundary layer Options:
max_target_frequency = 1e8
boundary_layer_material = material_name[0]
number_of_layers = 3



# Generating OCC primitive boxes:
bar1 = Box(Pnt(-1,0,0), Pnt(0,1,1))
bar2 = Box(Pnt(0,0,0), Pnt(1,1,1))

bar1 = bar1.Rotate(Axis((0,0,0),Y), 33)
bar2 = bar2.Rotate(Axis((0,0,0),Y), 33)


# Generating surrounding non-conducting region as [-1000,1000]^3 box:
outer_box = Box(Pnt(-1000, -1000, -1000), Pnt(1000,1000,1000))

# setting material and bc names:
# For compatability, we want the non-conducting region to have the 'outer' boundary condition and be labeled as 'air'
bar1.mat(material_name[0])
bar2.mat(material_name[1])
bar1.bc('default')
bar2.bc('default')
outer_box.mat('air')
outer_box.bc('outer')

# Setting maxh for each object:
bar1.maxh = 0.12
bar2.maxh = 0.12
outer_box.maxh = 100

# Joining the three meshes:
# Glue joins two OCC objects together without interior elements.
joined_object = Glue([bar1, bar2, outer_box])

# Generating Mesh:
geo = OCCGeometry(joined_object)
nmesh = geo.GenerateMesh()

# Applying Boundary Layers:
mu0 = 4 * 3.14159 * 1e-7
tau = (2/(max_target_frequency * sigma[1] * mu0 * mur[0]))**0.5 / alpha
layer_thicknesses = [(2**n)*tau for n in range(number_of_layers)]

nmesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=boundary_layer_material,
                           domains=boundary_layer_material, outside=False)


nmesh.Save(r'VolFiles/OCC_dualbar_rotated.vol')
ngmesh = Mesh(nmesh)

