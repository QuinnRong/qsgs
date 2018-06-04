# 2D/3D poisson solver for composite material
#
# First added:  Mostafa Mollaali 2018-04-27

from dolfin import *
import random
from os import listdir
from os.path import isfile, join
import numpy as np
import os
import time

def add_epsilon(data):
    shape = data.shape
    print(shape)
    num_zero = 0
    for index in range(data.size):
        if (data[index] < epsilon):
            data[index] = epsilon
            num_zero += 1
    print(num_zero)

def pre_process(data):
    print(data.shape)
    volume_fraction = 0
    for index in range(data.size):
        if (data[index] < 0.5):
            data[index] = TC[0]
        else:
            data[index] = TC[1]
            volume_fraction += 1

#---------------------------------------------------------------------------------
# Sub domain for Dirichlet boundary condition
#---------------------------------------------------------------------------------
class Top_z1(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[ndim-1] - 1.0) < DOLFIN_EPS and on_boundary

class Bottom_z0(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[ndim-1] -0.0) < DOLFIN_EPS and on_boundary

class Right_x1(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] - 1.0) < DOLFIN_EPS and on_boundary

class Left_x0(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] -0.0) < DOLFIN_EPS and on_boundary


imported_kappa = True
dim = 50
print(dim)
TC = [1, 10]
epsilon = 1e-12
ndim = 3 #2D or 3D problem

for idx in range(0, 100):
    dump_file = open("../structure/rand-50-0.250000/result.txt", "a")
#----------------------------------------------------------------------------
# import the kappa
#----------------------------------------------------------------------------
    time_start=time.time()
    if imported_kappa:
        dirname = "../structure/rand-50-0.250000/" + str(int(idx))
        #check directory exists
        if(os.path.exists(dirname)):
            print("Directory Exists")
            kappa_input = []
            for j in range(0, dim): 
                filename = dirname+"/3D_"+str(int(idx))+"_"+str(int(j))+".dat"
                #check file exists
                if(os.path.exists(filename)):
                    #print("File Exists")
                    # xx = np.loadtxt(filename,  unpack=True) # open(filename)
                    kappa=  np.loadtxt(filename,  unpack=True) .ravel()
                    kappa_input = np.append(kappa_input, kappa)
                    length=len(kappa_input)

                else:
                    print("File does not exists")
        else:
            print("Directory does not exists")
    if imported_kappa:
        pre_process(kappa_input)
#----------------------------------------------------------------------------
# Create mesh 
#----------------------------------------------------------------------------
    if ndim == 2:
            mesh = UnitSquareMesh.create(32,32, CellType.Type_quadrilateral)
    elif ndim == 3 :
        mesh = UnitCubeMesh.create(dim,dim,dim, CellType.Type_hexahedron)
    ndim = mesh.geometry().dim() # get number of space dimensions
#----------------------------------------------------------------------------
# Create mesh functions for conductivity (k00
#----------------------------------------------------------------------------
    k00 = MeshFunction("double", mesh, mesh.geometry().dim())

# Iterate over mesh and set values
    for (i,cell) in enumerate(cells(mesh)):

        
        if imported_kappa:
            k00[cell] = kappa_input[i] 

        else:
            k00[cell] = 1.# cell.midpoint().z()
            #if random.random() < 0.5:
            #   k00[cell] = 10.0

            #else:
            #   k00[cell] = 1.0


# Store to file
    mesh_file = File("UnitMesh.xml.gz")
    k00_file = File("UnitMesh_k00.xml.gz")
    k00_file << k00
    File("k00.pvd") << k00
#----------------------------------------------------------------------------
# create function space
#----------------------------------------------------------------------------
    V = FunctionSpace(mesh, "Lagrange", 1)
#----------------------------------------------------------------------------
# Code for C++ evaluation of conductivity
#----------------------------------------------------------------------------
    conductivity_code = """
    class Conductivity : public Expression
    {
    public:

      // Create expression with 3 components
      Conductivity() : Expression(1) {}

      // Function for evaluating expression on each cell
      void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
      {
        const uint D = cell.topological_dimension;
        const uint cell_index = cell.index;
        values[0] = (*k00)[cell_index];

      }

      // The data stored in mesh functions
      std::shared_ptr<MeshFunction<double> > k00;

    };
    """
    c = Expression(cppcode=conductivity_code, degree=0)
    c.k00 = k00

    kappa= c[0]

    # Initialize sub-domain instances
    top = Top_z1()
    bottom = Bottom_z0()
    right= Right_x1()
    left= Left_x0()

# Initialize mesh function for boundary domains
    boundaries = MeshFunction("size_t", mesh, mesh.geometry().dim()-1)
    boundaries.set_all(0)
    top.mark(boundaries, 1)
    bottom.mark(boundaries, 2)
    right.mark(boundaries, 3)
    left.mark(boundaries, 4)

# Define boundary condition
    bc_1 = DirichletBC(V, Constant(1.0), boundaries, 1)
    bc_2 = DirichletBC(V, Constant(0.0), boundaries, 2)
    bc=[bc_1, bc_2]
# Define measure for boundary condition integral
    ds = Measure("ds")(subdomain_data=boundaries)
#---------------------------------------------------------------------------------
# Define variational problem
#---------------------------------------------------------------------------------
    u = Function(V)
    v = TestFunction(V)

    f = Constant(0.0)
    F = dot(kappa*grad(u), grad(v))*dx - f*v*dx
#---------------------------------------------------------------------------------
# Compute solution
#---------------------------------------------------------------------------------
#solve(F == 0, u, bc,solver_parameters={"newton_solver":{"relative_tolerance":1e-6}})
#solve(F == 0, u, bc,solver_parameters={"newton_solver":{"relative_tolerance":1e-6},
#				       "newton_solver":{"absolute_tolerance":1e-8},
#				       "newton_solver":{"linear_solver":"mumps"}})
    solve(F == 0, u, bc,solver_parameters={"newton_solver":{"relative_tolerance":1e-4},
				       "newton_solver":{"absolute_tolerance":1e-6},
				       "newton_solver":{"linear_solver":"mumps"}})
#---------------------------------------------------------------------------------
# Post proccessing
#---------------------------------------------------------------------------------



# Evaluate integral of normal gradient over top boundary
    n = FacetNormal(mesh)

    m1 = kappa*dot(grad(u), n)*ds
    v1 = assemble(m1)
    print("kappa * grad(u) * n ds = ", v1)

    Flux_Top = kappa*dot(grad(u), n)*ds(1)
    Total_Flux_Top = assemble(Flux_Top)
    print("kappa * grad(u) * n ds(1) = ", Total_Flux_Top)

    Flux_bottom = kappa*dot(grad(u), n)*ds(2)
    Total_Flux_bottom = assemble(Flux_bottom)
    print("kappa * grad(u) * n ds(2) = ", Total_Flux_bottom)

    Flux_right = kappa*dot(grad(u), n)*ds(3)
    Total_Flux_right = assemble(Flux_right)
    print("kappa * grad(u) * n ds(3) = ", Total_Flux_right)


    Flux_left = kappa*dot(grad(u), n)*ds(4)
    Total_Flux_left = assemble(Flux_left)
    print("kappa * grad(u) * n ds(4) = ", Total_Flux_left)
    # Save solution in VTK format
    file = File("poisson.pvd")
    file << u

    dump_file.write("%4d%10.6f\n" % (idx, (Total_Flux_Top - Total_Flux_bottom) / 2))
    dump_file.close()

    time_end=time.time()
    print('totally cost',time_end-time_start)
