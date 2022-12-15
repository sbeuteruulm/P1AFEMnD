using SparseArrays
using LinearAlgebra
using Combinatorics
using GroupSlices
using PyPlot
using DelimitedFiles
using LinearOperators
using Preconditioners
using IterativeSolvers
using MAT
include("dO.jl")
include("volngrad.jl")
include("provideHyperFaceData.jl")
include("quadnd.jl")
include("solvePDE.jl")
include("provideNeighbourData.jl")
include("createSimultBisection.jl")
include("refine.jl")
include("refineNVBnD.jl")
include("computeEtaR.jl")
include("adaptiveAlgorithm.jl")
include("extractBoundary.jl")

# Define 2d mesh with L-shaped domain ([-1,1]x[-1,1])/([0,1]x[0,1])
coordinates = [-1.0 -1.0;
                0.0 -1.0;
               -1.0  0.0;
                0.0  0.0;
                1.0  0.0;
               -1.0  1.0;
                0.0  1.0;
                1.0  1.0]

elements = [1 2 4;
            1 3 4;
            6 3 4;
            6 7 4;
            4 7 8;
            4 5 8]

dirichlet = [2 1;
             1 3;
             3 6;
             6 7;
             7 8;
             8 5]

neumann = [2 4;
           4 5]

level = zeros(Int64,size(elements,1),1);

# Define problem specification
D = dO(1.0*Matrix(I,2,2),zeros(2,1),0)

# Define right hand side
function f(x)
	return ones(size(x,1),1)
end

function uD(x)
	return zeros(size(x,1),1)
end

function g(x)
	return zeros(size(x,1),1)
end

# Define parameters for stopping criterion
nEmax = 10000
eta = 0.7

# Run FEM
x,energy,coordinates,elements,bdry = adaptiveAlgorithm(coordinates,elements,
    level,f,g,uD,D,nEmax,eta,dirichlet,neumann)

# Plot solution
plot_trisurf(coordinates[:,1],coordinates[:,2], x[:],
    triangles=elements-ones(size(elements)),cmap=plt.cm.jet,linewidth=100);
