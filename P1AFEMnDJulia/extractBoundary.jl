function extractBoundary(elements,element2neighbour,nB)
    # Extract boundary data from a FEM mesh with given information of elements
    # and neighbour relations
    #
    # boundary = extractBoundary(elements,element2neighbour,nB) boundary is
    # a tuple with nB boundary arrays containing the boundary hyperfaces
    # by their vertex numbers row wise. The information is extracted from the
    # array elements that saves the vertex numbers of each simplex and the
    # array element2neighbour that contains the element numbers of all
    # neighbours per element. A negative number in element2neighbour implies a
    # boundary type.
    #
    # Comments:
    #   During the refinement process element and neighbour information are
    #   updated but not the boundary arrays. They can be identified by the
    #   entries in the element2neighbour array. And thus the refined boundary
    #   is given and can be used for the approximation of a solution.
    #
    # Remark:
    #   This program is a supplement to the paper
    #   >> Efficient P1-FEM for any space dimension in Matlab <<
    #   by S. Beuter, and S. Funken. The reader should
    #   consult that paper for more information.
    #
    # Authors:
    #   S. Beuter, S. Funken 18-10-22

    nE = size(elements,1)
    nN = size(elements,2)
    facets = vcat(collect(combinations(1:nN,nN-1))[end:-1:1]'...) .- 1
    boundary = ()
    for k = 1:nB
        rowCol = findall(element2neighbour .== -k)
        idx = getindex.(rowCol,1) * ones(Int64,1,nN-1) +
            facets[getindex.(rowCol,2),:] .* nE
        if !isempty(rowCol)
            boundary = tuple(boundary...,(reshape(elements[idx],
                length(rowCol),:),)...)
        else
            boundary = tuple(boundary...,(zeros(0,nN-1),)...)
        end
    end
    return boundary
end
