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
include("nFicheraCube.jl")
include("nSubCube.jl")
include("KuhnTriangulation.jl")
include("convertdec.jl")

# Define settings
    # problem size ("small", "medium", "large")
    Psize = "medium"
    # adaptivity ("low", "high", "uniform")
    adaptivity = "high"
    # problem dimesnision (a natural number)
    d = 3
    # problem type ("linear", "quadratic" (in the first component))
    pType = "quadratic"

# Define stopping criterion and adaptivity parameter according to problem
# setting
if Psize == "small" && adaptivity == "low"
    nEmax = 1000
    eta = 0.7
elseif Psize == "small" && adaptivity == "high"
    nEmax = 1000
    eta = 0.2
elseif Psize == "small" && adaptivity == "uniform"
    nEmax = 1000
    eta = 1
elseif Psize == "medium" && adaptivity == "low"
    nEmax = 50000
    eta = 0.7
elseif Psize == "medium" && adaptivity == "high"
    nEmax = 50000
    eta = 0.2
elseif Psize == "medium" && adaptivity == "uniform"
    nEmax = 50000
    eta = 1
elseif Psize == "large" && adaptivity == "low"
    nEmax = 1000000
    eta = 0.7
elseif Psize == "large" && adaptivity == "high"
    nEmax = 1000000
    eta = 0.2
elseif Psize == "large" && adaptivity == "uniform"
    nEmax = 1000000
    eta = 1
else
    println("No possible setting choice")
end

# Define mesh dependent on dimension d and subcube size 1/(s+1). The domain
# is a Fichera cube [-1,1]^d/[0,1]^d. Here the whole boundary is formed of
# dirichlet boundary, while the neumann boundary is empty.
s = 1
coordinates,elements,dirichlet = nFicheraCube(d,s,"tucker")
neumann = []

level = zeros(Int64,size(elements,1),1)

# Define problem specification
D = dO(d.*Matrix(I,d,d)+ones(d,d),zeros(d,1),0)

# Define right hand side
if pType == "linear"
    function f(x)
    	return zeros(size(x,1),1)
    end
    function g(x)
    	return zeros(size(x,1),1)
    end
    function uD(x)
    	return x[:,1]
    end
elseif pType == "quadratic"
    function f(x)
    	return -(2*d+2).*ones(size(x,1),1)
    end
    function g(x)
    	return zeros(size(x,1),1)
    end
    function uD(x)
    	return x[:,1].^2 + x[:,2]
    end
else
    disp("problem type not supported.")
end

# Run FEM
x,energy,coordinates,elements,bdry = adaptiveAlgorithm(coordinates,elements,
    level,f,g,uD,D,nEmax,eta,dirichlet,neumann)

# Relative error
if (pType == "linear")
    maxi,idx = findmax(abs.(x[:]-coordinates[:,1]))
    println(maxi)
    if coordinates[idx,1] != 0
        println("Maximum relative error ", string(abs(maxi)/
            abs(coordinates[idx,1])))
    end
else
    maxi,idx = findmax(abs.(x[:]-(coordinates[:,1].^2+coordinates[:,2])))
    if (coordinates[idx,1].^2+coordinates[idx,2]) != 0
        println("Maximum relative error ", string(abs(maxi)/
            abs(coordinates[idx,1].^2+coordinates[idx,2])))
    end
end

if size(coordinates,1) < 1e5
    if d == 2
        plot_trisurf(coordinates[:,1],coordinates[:,2], x[:],
            triangles=elements-ones(size(elements)),cmap=plt.cm.jet,linewidth=100)
    end
end
