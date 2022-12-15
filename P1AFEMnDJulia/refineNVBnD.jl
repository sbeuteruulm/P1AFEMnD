function refineNVBnD(coordinates,elements,element2neighbour,level,marked,boundary)
    # This function gives a framework for the refinement process.
    #
    # coordinates,elements,element2neighbour,level = ...
    # refineNVBnD(coordinates,elements,element2neighbour,level,marked,boundary)
    # updates a mesh by the newest vertex bisection. The mesh is described by
    # the arrays coordinates containing the coordinates for each vertex row
    # wise, the array elements that describes the simplices by their vertex
    # numbers, element2neighbour assembles for each element all its neighbours
    # in one row, and vector level that gives the level of each simplex. These
    # arrays get enlarged within the function before they are returned. Vector
    # marked contains the numbers off all elements that are marked for
    # refinement. The tuple boundary assembles all boundary arrays.
    #
    # Comments:
    #   This function is called to perform the REFINE step of the adaptive FEM.
    #   All marked elements are bisected by the newest vertex bisection. To
    #   ensure that every marked element and neighbours sharing the refinement
    #   edge are bisected exactly once, the refine function is called if the
    #   global level has not been updated before. The refine function
    #   identifies all neighbours sharing a refinement edge and initiates the
    #   Newest Vertex Bisection.
    #
    # Remark:
    #   This program is a supplement to the paper
    #   >> Efficient P1-FEM for any space dimension in Matlab <<
    #   by S. Beuter, and S. Funken. The reader should
    #   consult that paper for more information.
    #
    # Authors:
    #   S. Beuter, S. Funken 18-10-22

    global C, E, L, E2N, i2fi, nC, nE
    nB = length(boundary)
    E2N = element2neighbour
    i2fi = 0
    nE = size(elements,1)
    nC = size(coordinates,1)
    C = coordinates
    E = elements
    L = level
    # refine mesh by bisection
    lev = L
    for e = marked
        if lev[e] == L[e]
            refine(e)
        end
    end
    #return result
    element2neighbour = E2N
    E2N = E2N[1:nE,:]
    return C[1:nC,:], E[1:nE,:], element2neighbour, L[1:nE]
end
