function provideHyperFaceData(elements,bdry)
    # Calculation of additional mesh data. Especially, the hyperface data is
    # collected.
    #
    # hyperface2nodes,element2hyperfaces,boundary2hyperfaces = ...
    # provideHyperFaceData(elements,bdry) Creation of a numbering for
    # hyperfaces of a mesh described by their node numbers, which are given for
    # each element in elements. Each hyperface is given by its nodes row wise
    # in hyerface2elements. The hyperfaces of the simplices are assembled in
    # element2hyperfaces where each row contains the hyperfaces of an element.
    # According to the tuple bdry which contains all nodes of each
    # boundary type, a cell array boundary2hyperfaces is created that contains
    # all hyperfaces of each boundary type.
    #
    # Comments:
    #   For any boundary type the hyperfaces being part of it are also given in
    #   an array.
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
    nB = length(bdry)
    nbE = zeros(Int64,nB)
    for j = 1:nB
        nbE[j] = size(bdry[j],1)
    end
    ptr = cumsum([nE*nN; nbE])
    # i'th hyperface is opposite the i'th node
    ordering = vcat(collect(combinations(1:nN,nN-1))[end:-1:1]'...)
    # ascending sorting of nodes of all facets
    facets = [sort(reshape(elements[:,ordering],:,nN-1),dims=2)
        zeros(Int64,ptr[end]-ptr[1],nN-1)]
    # add boundary facets
    for j = 1:nB
        if !isempty(bdry[j])
            facets[ptr[j]+1:ptr[j+1],:] = sort(bdry[j],dims=2)
        end
    end
    # create unique face numbering, element2hyperfaces, and boundary2hyperfaces
    hyperface2nodes = unique(facets,dims=1)
    uniqueIdx = groupinds(groupslices(facets,dims=1))
    sortIdx = sortslices([hyperface2nodes 1:size(hyperface2nodes,1)],dims=1)
    hyperface2nodes = sortIdx[:,1:end-1]
    sortIdx = sortslices([sortIdx[:,end] 1:size(sortIdx,1)],dims=1)[:,end]
    J = zeros(Int64,size(facets,1))
    for k = 1:length(uniqueIdx)
       J[uniqueIdx[k]] .= sortIdx[k]
    end
    element2hyperfaces = reshape(J[1:nN*nE],:,nN)
    result = (hyperface2nodes, element2hyperfaces)
    for j = 1:nB
        result = tuple(result...,
            (reshape(J[ptr[j]+1:ptr[j+1]],nbE[j],1),)...)
    end
    return result
end
