from netgen.occ import *
# from ngsolve import *
from netgen.meshing import BoundaryLayerParameters
import math

"""
Paul Ledger - 2025:
coil example for Netgen OCC geometry mesh generation.
Object has prismatic boundary layer elements added.

EDIT 2023:
Netgen-Mesher version 6.2.2301 gives a different result for the assigned materials when compared to version 6.2.2204.
The material assinged to 'box' should be 'air', and indeed this is what is reported when using the older version of
netgen. When using the new version, it reports the material as 'default'.

To test this I uninstalled both ngsolve and netgen-mesher and reinstalled both using the command
pip3 install ngsolve==6.2.2204

"""



# Setting mur, sigma, alpha, and defining the top level object name:
material_name = ['mat1']
mur = [100]
sigma = [6e6]
alpha = 0.01

# Boundary Layer Settings: max frequency under consideration, the total number of prismatic layers and the material of each layer.
# Setting Boundary layer Options:
max_target_frequency = 1e8
boundary_layer_material = material_name[0]
number_of_layers = 2

# radius of the wire
rcoil=0.05
# radius and height of the helix
rheli=0.25
hheli=2
cyl = Cylinder((0,0,0), Z, r=rheli, h=hheli).faces[0]
heli = Edge(Segment((0,0), (20*math.pi, hheli)), cyl)
ps = heli.start
vs = heli.start_tangent
pe = heli.end
ve = heli.end_tangent


# create a wire model and then pipe a wire of radius rcoil along its length
spiral = Wire([heli])
circ = Face(Wire([Circle(ps, vs, rcoil)]))
coil = Pipe(spiral, circ)

coil.faces.maxh=10
coil.mat(material_name[0])
coil.bc("default")

# Generating a large non-conducting region. For compatability with MPT-Calculator, we set the boundary condition to 'outer'
# and the material name to 'air'.
box = Box(Pnt(-100, -100, -100), Pnt(100,100,100))
box.bc('outer')
box.maxh=100
box.mat("air")
air=box-coil
air.mat('air')
# Joining the two meshes:
# Glue joins two OCC objects together without interior elemements
joined_object = Glue([coil, air])
print("got geometry")


# Generating Mesh:
#nmesh = OCCGeometry(joined_object).GenerateMesh(meshsize.coarse, maxh=100)

print("got mesh")

# Creating Boundary Layer Structure:
if number_of_layers > 0:
    mu0 = 4 * 3.14159 * 1e-7
    tau = (2/(max_target_frequency * sigma[0] * mu0 * mur[0]))**0.5 / alpha
    layer_thicknesses = [(2**n)*tau for n in range(number_of_layers)]

    B = BoundaryLayerParameters(boundary=".*", thickness=layer_thicknesses, new_material=boundary_layer_material,
                           domain=boundary_layer_material, outside=False)#, disable_curving=False)

    nmesh = OCCGeometry(joined_object).GenerateMesh(meshsize.coarse, boundary_layers=[B])
    #nmesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=boundary_layer_material,
    #                       domains=boundary_layer_material, outside=False)


nmesh.Save(r'VolFiles/OCC_coil_al_0_01_r_0_25_h_2.vol')
# print(nmesh.GetMaterial(2))
from ngsolve import *
mesh = Mesh(nmesh)
print(f'Materials = {mesh.GetMaterials()}')
