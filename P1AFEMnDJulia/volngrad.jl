function volngrad(coordinates,elements,flag)
    # Calculation of volumes of each element and gradients of nodal interpolant
    # functions
    #
    # vol, grad = volngrad(coordinates,elements,flag) or
    # vol = volngrad(coordinates,elements,flag) calculates the volumes and
    # gradient of the basis functions for a given simplex. If flag is set 2 then
    # two output parameters vol and grad are expected and returned. The latter
    # is optional and disregareded if flag is any other value than 2. The
    # gradients are given as a Matrix of arrays for each element. The input
    # array coordinates saves for any mesh node the respective coordinates in a
    # row. In the same way, the nodes of each element are given in the
    # respective row of elements.
    #
    # Comments:
    #   A Cholesky decomposition of the matrices containing the differences of
    #   the element vertices is calculated. The advantage is that the volumes
    #   can be calculated as a byproduct.
    #
    # Remark:
    #   This program is a supplement to the paper
    #   >> Efficient P1-FEM for any space dimension in Matlab <<
    #   by S. Beuter, and S. Funken. The reader should
    #   consult that paper for more information.
    #
    # Authors:
    #   S. Beuter, S. Funken 18-10-22

    nD = size(coordinates,2)
    nN = size(elements,2)
    nE = size(elements,1)
    if nN == 2
        volume = abs.(coordinates[elements[:,1]]-coordinates[elements[:,2]])
        if flag == 2
            grad = Matrix(undef,nN,1)
            grad[1] = 1 ./(coordinates[elements[:,1]]-coordinates[elements[:,2]])
            grad[2] = 1 ./(coordinates[elements[:,2]]-coordinates[elements[:,1]])
            return volume, grad
        end
        return volume
    else
        # Compute Gram matrix L
        idx = repeat(1:nN-1,1,nN-1)
        ii = nonzeros(sparse(UpperTriangular(idx)'))
        jj = nonzeros(sparse(LowerTriangular(idx)))
        grad = Matrix(undef,nN,1)
        grad[nN] = coordinates[elements[:,nN],:]
        for k = 1:nN-1
            grad[k] = coordinates[elements[:,k],:] - grad[nN]
        end
        L = zeros(nE,nN*(nN-1)÷2)
        for k = 1:length(ii)
            L[:,k] .+= sum(grad[ii[k]].*grad[jj[k]],dims=2)[:]
        end
        # Compute Cholesky decomposition of L
        ptr = [1; 1 .+ cumsum(nN-1:-1:1)]
        L[:,ptr[1]] = sqrt.(L[:,ptr[1]])
        for k = 1:nN-2
            L[:,ptr[k]+1:ptr[k+1]-1] ./= L[:,ptr[k].*ones(Int64,
                ptr[k+1]-ptr[k]-1)]
            for j = k+1:nN-1
                L[:,ptr[j]:ptr[j+1]-1] .-= L[:,ptr[k]+j-k:ptr[k+1]-1] .*
                    L[:,ones(Int64,nN-j).*(ptr[k].+j.-k)]
            end
            L[:,ptr[k+1]] = sqrt.(L[:,ptr[k+1]])
        end
        # Compute volume
        volume = prod(L[:,ptr[1:(nN-1)]],dims=2)./factorial(nN-1)
        if flag == 2
            # Compute L y = A'
            for j = 1:nN-2
                grad[j] ./= L[:,ptr[j].*ones(Int64,nD)]
                for k = j+1:nN-1
                    grad[k] .-= grad[j].*L[:,(ptr[j].+k.-j).*ones(Int64,nD)]
                end
            end
            grad[nN-1] ./= L[:,ptr[nN-1].*ones(Int64,nD)]
            #  Compute L' * x = y;
            grad[nN-1] ./= L[:,ptr[nN-1].*ones(Int64,nD)]
            grad[nN] = -grad[nN-1]
            for j = nN-2:-1:1
                for k = j+1:nN-1
                    grad[j] .-= grad[k] .* L[:,(ptr[j].+k.-j).*ones(Int64,nD)]
                end
                grad[j] ./= L[:,ptr[j.*ones(Int64,nD)]]
                grad[nN] .-= grad[j]
            end
            return volume, grad
        end
    end
    return volume
end
