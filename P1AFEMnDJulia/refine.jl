function refine(e)
    # Performance of closure step of Newest Vertex Bisection refinement.
    #
    # s = refine(e) for a given element number e that shall be refined,
    # additional neighbouring elements have to be refined, to keep the mesh
    # conforming. All numbers of elements that get additionally refined are
    # returned in s.
    #
    # Comments:
    #   Necessary mesh information is given as global arrays. For a marked
    #   element any neighbour sharing the refinement edge has to be refined
    #   as well. If two elements do not have the refinement edge in common,
    #   additional refinements may be necessary. Then, bisection of all
    #   elements around a refinement edge is done by the call of
    #   createSimultBisection.
    #
    # Remark:
    #   This program is a supplement to the paper
    #   >> Efficient P1-FEM for any space dimension in Matlab <<
    #   by S. Beuter, and S. Funken. The reader should
    #   consult that paper for more information.
    #
    # Authors:
    #   S. Beuter, S. Funken 18-10-22

    global E, E2N, nE
    nN = size(E,2)
    nD = nN-1
    K = zeros(Int64,0,1)
    F = e
    FK = zeros(Bool,nE)
    FK[e] = true
    while !isempty(F)
        Fnew = zeros(Int64,0,1)
        for e1 = F
            e1a = E[e1,1]
            e1b = E[e1,nN]
            E2NLocal = E2N[e1,1:end]
            for j = 2:nD
                e2 = E2NLocal[j]
                if e2 > 0 && !FK[e2]
                    # e2 has same refinement edge as e1
                    e2a = E[e2,1]
                    e2b = E[e2,nN]
                    if ((e1a == e2a) && (e1b==e2b)) ||
                        ((e1a == e2b) && (e1b == e2a))
                        Fnew = [Fnew;e2]
                    else
                        s = refine(e2)
                        FK = [FK;zeros(Bool,nE-length(FK))]
                        # add to Fnew the child of e2 that is a neighbour of e1
                        if sum(in.(s[1],E2N[e1,1:end])) > 0
                            Fnew = [Fnew;s[1]]
                        else
                            Fnew = [Fnew;s[2]]
                        end
                    end
                end
            end
        end
        K = [K;F]
        F = Fnew
        FK[F] .= true
    end
    s = [K[1] nE+1]
    K = sort(unique(K))
    i = findall(K.==e)
    K = [K[i[1]:end]; K[1:(i[1]-1)]]
    createSimultBisection(K)
    return s
end
