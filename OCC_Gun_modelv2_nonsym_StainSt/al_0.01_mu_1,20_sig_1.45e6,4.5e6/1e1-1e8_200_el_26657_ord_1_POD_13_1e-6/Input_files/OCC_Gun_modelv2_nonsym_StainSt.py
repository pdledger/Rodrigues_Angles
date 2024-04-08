
from netgen.occ import *
from netgen.webgui import Draw as DrawGeo

# Set materials
material_name = ['Steel','StainSteel']
sigma = [4.5e6, 1.45e6]
mur = [20, 1]# reduced mU-r for testing was [100,1]
alpha =0.01

# Boundary layer options
max_target_frequency = 1e8
number_of_layers = 3 

# Reciever
box = Box(Pnt(12,-0.5,-15.75), Pnt(20,0.5,-0.75))
# Barrel
cyl = Cylinder(Pnt(0,0,0), X, r=2, h=20)
cyl_inside = Cylinder(Pnt(0,0,0), X, r=1, h=18)
# Barrel-extra
barrel_ex = Box(Pnt(0,-1.5,-0.5),(20,-2.5,0.5))


# Reciever is the box without the barrel
barrel_ex_no_cyl = barrel_ex - cyl
barrel=cyl-cyl_inside + barrel_ex_no_cyl
rec=box-cyl
# Set materials
barrel.mat(material_name[0])
rec.mat(material_name[1])

# Set bcs
barrel.bc("gunbc")
rec.bc("gunbc")

barrel.maxh = 1
rec.maxh = 1
# Join objects together to form the gun
gun = Glue([rec,barrel])
DrawGeo (gun)

# Add large outer box
box = Box(Pnt(-1000, -1000, -1000), Pnt(1000,1000,1000))
box.mat('air')
box.bc('outer')
box.maxh = 1000

# Join gun to the box
joined_object = Glue([gun, box])
geo = OCCGeometry(joined_object)
nmesh = geo.GenerateMesh(grading=0.5)

# Setup the boundary layers
mu0 = 4 * 3.14159 * 1e-7
s = alpha / 1e-3
tau = (2/(max_target_frequency * sigma[0] * mu0 * mur[0]))**0.5 /(alpha/s)
print(s,tau)
layer_thicknesses = [ (2**n)*tau for n in range(number_of_layers)]

nmesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=material_name[0],
                           domains=material_name[0], outside=False)

#tau = (2/(max_target_frequency * sigma[1] * mu0 * mur[1]))**0.5 /alpha
#layer_thicknesses = [ (2**n)*tau for n in range(number_of_layers)]

#nmesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=material_name[1],
#                           domains=material_name[1], outside=False)

#Save the mesh
print("save mesh")
nmesh.Save(r'VolFiles/OCC_Gun_modelv2_nonsym_StainSt.vol')

