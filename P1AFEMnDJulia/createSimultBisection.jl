function createSimultBisection(marked)
    # Simultaneous bisection of a marked element and its neighbours sharing
    # the refinement edge.
    #
    # createSimultBisection(marked) updates a globally defined mesh by
    # bisection of a marked element and its neighbours sharing a refinement
    # edge. An element number of the simplex that shall be refined is expected
    # as input marked.
    #
    # Comments:
    #   All elements around the refinement edge can be bisected
    #   simultaneously. The element array gets enlarged. Since the
    #   neighbour information is crucial for the refinement process, it has to
    #   be updated for all bisected elements, too. With each call of
    #   createSimultBisection one new mesh node is calculated and added to the
    #   array coordinates.
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
    nR = length(marked) # number of simplices to be refined
    nN = size(E,2)      # number of nodes per simplex
    # create new node
    ae = E[marked[1],[1,end]]
    C = [C; sum(C[ae,:],dims=1)/2]
    nC += 1
    # create new elements and update element2neighbour
    elem = E[marked,:]
    low = [elem[:,1] ones(Int64,nR,1).*nC elem[:,2:end-1]]
    up = [elem[:,end] ones(Int64,nR,1).*nC elem[:,2:end-1]]
    neigh = E2N[marked[:],1:end]
    n_low = [nE.+(1:nR) neigh[:,end] neigh[:,2:end-1]]
    i2fi = [i2fi;zeros(maximum([0,maximum(marked)-length(i2fi)]))]
    i2fi[marked] = 1:nR
    tmp = neigh[:,2:end-1]
    idx = findall(tmp[:].>0)
    tmp[idx] = nE.+i2fi[tmp[idx]]
    n_up = [marked[:] neigh[:,1] tmp]
    ind = elem[:,1] .> elem[:,end]
    lev = L[marked[:]].%(nN-1)
    tmp = n_up[ind,2:end]
    n_up[ind,2:end] = n_low[ind,2:end]
    n_low[ind,2:end] = tmp
    # reestablish correct node order
    for k = 0:nN-3 #nD-1
        idx = lev .== k
        up[idx,3+k:end] = up[idx,end:-1:3+k]
        n_up[idx,3+k:end] = n_up[idx,end:-1:3+k]
    end
    # update global arrays for elements and nieghours
    tmp = n_up[ind,3:end]
    n_up[ind,3:end] = n_low[ind,3:end]
    n_low[ind,3:end] = tmp
    E2N[marked,:] = n_low
    E2N = [E2N; n_up]
    tmp = up[ind,:]
    up[ind,:] = low[ind,:]
    low[ind,:] = tmp
    E[marked,:] = low
    E = [E; up]
    # update element2neighbour for elements opposite the first node
    ii = findall(neigh[:,1] .> 0)
    imat = E2N[neigh[ii,1],:] .==
        marked[ii] * ones(Int64,1,nN)
    posIndex = findall(imat.>0)
    tmp = n_low[:,1]
    tmp[ind] = n_up[ind,1]
    for k = 1:length(posIndex)
        E2N[neigh[ii[posIndex[k][1]],1],posIndex[k][2]] =
            tmp[ii[posIndex[k][1]]]
    end
    # update element2neighbour for elements opposite the last node
    ii = findall(neigh[:,end] .> 0)
    imat = E2N[neigh[ii,end],:] .==
        marked[ii] * ones(Int64,1,nN)
    posIndex = findall(imat.>0)
    tmp = n_up[:,1]
    tmp[ind] = n_low[ind,1]
    for k=1:length(posIndex)
        E2N[neigh[ii[posIndex[k][1]],end],posIndex[k][2]] =
            tmp[ii[posIndex[k][1]]]
    end
    # update level
    L = [L; zeros(Int64,nR)]
    L[marked] .+= 1
    L[nE.+(1:nR)] .= L[marked]
    nE = nE + nR
end
