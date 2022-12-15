function nFicheraCube(n,s,pattern)
    # Calculation of mesh data for an n-dimensional Fichera cube, i.e. a
    # simplicial triangulation of [-1,1]^d / ([0,1]^d)
    #
    # coordinates,elements,boundary = nFicheraCube(n,s,pattern) creation of
    # mesh data for an n-Fichera cube of dimension n, where each subcube has an
    # edge length of 1/(s+1) and is split into Kuhn simplices. Facing subcubes
    # can either be reflected in any special direction if the parameter pattern
    # is 'tucker' or they are translated. This can be determined by
    # 'freudenthal'. The resulting mesh is given by the arrays coordinates
    # containing the coordinates of each node in its rows, elements that gives
    # all simplices by its vertices, and boundary that contains all boundary
    # facets by its nodes.
    #
    # Remark:
    #   This program is a supplement to the paper
    #   >> Efficient P1-FEM for any space dimension in Matlab <<
    #   by S. Beuter, and S. Funken. The reader should
    #   consult that paper for more information.
    #
    # Authors:
    #   S. Beuter, S. Funken 18-10-22

    coordinates, elements = nSubCube(n,2*s+1,pattern)
    coordinates = -1 .+ 2 .* coordinates
    # Delete simplices with midpoints > 0 (i.e. each component is positive)
    midpts = zeros(size(elements,1),n)
    for j = 1:n+1
        midpts .+= coordinates[elements[:,j],:]
    end
    midpts ./= n+1
    idx = .!iszero.(sum(midpts.<0,dims=2))
    elements = elements[idx[:],:]
    nodes = sort(unique(elements))
    i2fi = zeros(Int64,size(coordinates,1),1)
    i2fi[nodes] = 1:length(nodes)
    coordinates = coordinates[nodes,:]
    elements = i2fi[elements]
    # Extract boundary
    boundary = zeros(Int64,0,n)
    for k = 1:n+1
        bdry = elements[:,[1:k-1;k+1:n+1]]
        midpts = zeros(size(bdry,1),n)
        for j = 1:n
            midpts .+= coordinates[bdry[:,j],:]
        end
        midpts ./= n
        idx = (maximum(abs.(midpts),dims=2) .>= 1-1e-4) .|
            (minimum(midpts,dims=2) .>= -1e-4)
        boundary = [boundary; bdry[idx[:],:]]
    end
    return (coordinates,elements,boundary)
end
