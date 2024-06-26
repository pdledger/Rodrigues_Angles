import subprocess
import runpy
from netgen.csg import *

#Function definition which creates a mesh for a given .geo file
#Inputs -name of a .geo file (string)
#       -How fine the mesh should be (integer 1<=n<=5)
def Meshmaker(Geometry,Mesh):
    """
    Function to generate Netgen Mesh from .geo file and saves it as a similarly named .vol file

    :param Geometry: str path to the .geo file to be meshed.
    :param Mesh: int mesh granularity:
                 1 = very coarse
                 2 = coarse
                 1 = moderate
                 1 = fine
                 1 = very fine

    edit: James Elgy - 11 Oct 2022:
    Currently pip installations of netgen do not allow all command line arguments.
    See https://ngsolve.org/forum/ngspy-forum/1595-loading-geometry-from-command-line#4357

    I've added an option to mesh the .geo file using the CSGeometry package.

    """


    #Remove the .geo part of the file extention
    objname=Geometry[:-4]

    use_CSG = True
    if use_CSG:
        geo = CSGeometry('GeoFiles/' + Geometry)
        if Mesh == 'verycoarse':
            mesh = geo.GenerateMesh(meshsize.very_coarse)
        elif Mesh == 'coarse':
            mesh = geo.GenerateMesh(meshsize.coarse)
        elif Mesh == 'moderate':
            mesh = geo.GenerateMesh(meshsize.moderate)
        elif Mesh == 'fine':
            mesh = geo.GenerateMesh(meshsize.fine)
        elif Mesh == 'veryfine':
            mesh = geo.GenerateMesh(meshsize.very_fine)
        else:
            mesh = geo.GenerateMesh(maxh=Mesh)

        
        # material_name = ['tetra']
        # sigma = [1e6]
        # mur = [1]
        # alpha = 0.01

        # # Setting Boundary layer Options:
        # max_target_frequency = 1e8
        # boundary_layer_material = material_name[0]
        # number_of_layers = 1
        
        # mu0 = 4 * 3.14159 * 1e-7
        # tau = (2/(max_target_frequency * sigma[0] * mu0 * mur[0]))**0.5 / alpha
        # layer_thicknesses = [(2**n)*tau for n in range(number_of_layers)]

        # mesh.BoundaryLayer(boundary=".*", thickness=layer_thicknesses, material=boundary_layer_material, domains=boundary_layer_material, outside=False)
        
        mesh.Save('VolFiles/' + objname + '.vol')
        return

    #Define how fine the mesh will be
    if Mesh==1:
        Meshsizing='-verycoarse'
    elif Mesh==2:
        Meshsizing='-coarse'
    elif Mesh==3:
        Meshsizing='-moderate'
    elif Mesh==4:
        Meshsizing='-fine'
    elif Mesh==5:
        Meshsizing='-veryfine'
    else:
        print("No mesh created, please specify a number between 1-5")
    #Create the mesh
    try:
        subprocess.call(['netgen','-geofile=GeoFiles/'+Geometry,'-meshfile=VolFiles/'+objname+'.vol',Meshsizing,'-meshsizefile='+objname+'.msz','-batchmode'])
    except:
        subprocess.call(['netgen','-geofile=GeoFiles/'+Geometry,'-meshfile=VolFiles/'+objname+'.vol',Meshsizing,'-batchmode'])
    return



