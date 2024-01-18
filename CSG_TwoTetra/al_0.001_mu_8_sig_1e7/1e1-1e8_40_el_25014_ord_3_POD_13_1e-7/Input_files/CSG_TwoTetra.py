from netgen.csg import *


from netgen.occ import *

"""
James Elgy - 2023:
Two joined tetrahedra for Netgen OCC geometry loaded from existing step geometry.
"""

material_name = ['default']
sigma = [1e7]
mur = [8]
alpha = 0.001

# Setting Boundary layer Options:
max_target_frequency = 1e8
boundary_layer_material = material_name[0]
number_of_layers = 3


geo = CSGeometry(r'GeoFiles/TwoTetra.geo')

nmesh = geo.GenerateMesh(meshsize.coarse, optsteps3d=5, grading=0.6)
nmesh.SetMaterial(1, 'air')
nmesh.SetMaterial(2, 'default')
nmesh.SetMaterial(3, 'default')

nmesh.SetBCName(8, 'outer')




# Applying Boundary Layers:
mu0 = 4 * 3.14159 * 1e-7
tau = (2/(max_target_frequency * sigma[0] * mu0 * mur[0]))**0.5 / alpha
layer_thicknesses = [(2**n)*tau for n in range(number_of_layers)]

nmesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=boundary_layer_material, domains=boundary_layer_material, outside=False)

#print(layer_thicknesses)

nmesh.Save(r'VolFiles/CSG_TwoTetra.vol')