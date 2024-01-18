from netgen.occ import *

"""
James Elgy - 2023:
Two joined tetrahedra for Netgen OCC geometry loaded from existing step geometry.
"""

material_name = ['hex']
sigma = [1e7]
mur = [32]
alpha = 0.01

# Setting Boundary layer Options:
max_target_frequency = 1e8
boundary_layer_material = material_name[0]
number_of_layers = 3


geo = OCCGeometry(r'StepFiles/Irregular_hexahedron_concave.step')
hex = geo.shape.Move((-geo.shape.center.x, -geo.shape.center.y, -geo.shape.center.z))

print(hex.mass)

hex.bc('default')
hex.mat(material_name[0])
hex.maxh = 0.08

sph = Sphere(Pnt(0,0,0), r=100)
sph.mat('air')
sph.bc('outer')
sph.maxh=1000

joined_object = Glue([sph, hex])
nmesh = OCCGeometry(joined_object).GenerateMesh(meshsize.coarse)

# # Applying Boundary Layers:
mu0 = 4 * 3.14159 * 1e-7
tau = (2/(max_target_frequency * sigma[0] * mu0 * mur[0]))**0.5 / alpha
layer_thicknesses = [(2**n)*tau for n in range(number_of_layers)]

nmesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=boundary_layer_material, domains=boundary_layer_material, outside=False)


nmesh.Save(r'VolFiles/OCC_step_irreg_hexahedron_concave.vol')