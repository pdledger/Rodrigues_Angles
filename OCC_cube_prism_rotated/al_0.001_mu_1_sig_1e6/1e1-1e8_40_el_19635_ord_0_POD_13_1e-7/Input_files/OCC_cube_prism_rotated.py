from netgen.occ import *

# Setting mur, sigma, and defining the top level object name:
material_name = ['cube']
sigma = [1e6]
mur = [1]
alpha = 0.001

# Setting Boundary layer Options:
max_target_frequency = 1e8
boundary_layer_material = material_name[0]
number_of_layers = 3

cube = Box(Pnt(-5, -5, -5), Pnt(10, 5, 5))
# cube = cube.Rotate(Axis((0,0,0), Y), 33)

# Setting boundary conditions and material names.
cube.bc('default')
cube.mat(material_name[0])
cube.maxh = 1.5

# Generating a large non-conducting region. For compatability with MPT-Calculator, we set the boundary condition to 'outer'
# and the material name to 'air'.
box = Box(Pnt(-1000, -1000, -1000), Pnt(1000,1000,1000))
box.mat('air')
box.bc('outer')
box.maxh=1000

# Here we are joining the two geometries and generating the mesh.
joined_object = Glue([box, cube])
nmesh = OCCGeometry(joined_object).GenerateMesh()

# Applying Boundary Layers:
mu0 = 4 * 3.14159 * 1e-7
tau = (2/(max_target_frequency * sigma[0] * mu0 * mur[0]))**0.5 / alpha
layer_thicknesses = [(2**n)*tau for n in range(number_of_layers)]

nmesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=boundary_layer_material,
                           domains=boundary_layer_material, outside=False)

nmesh.Save(r'VolFiles/OCC_cube_prism_block.vol')
