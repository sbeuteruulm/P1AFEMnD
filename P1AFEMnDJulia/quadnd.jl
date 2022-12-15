function quadnd(nD,order)
    # Gauss type quadrature rules on n-dimensional simplex
    #
    # [lambda, weights] = quadnd(nD,order) gives the barycentric coordinates
    # and the weights for a quadrature formula of order order. The sum of
    # the weights is 1. The quadrature formula is exact for all polynomials
    # of total degree less or equalorder
    # (total degree = maximum of the degrees of all terms in a polynomial)
    #
    # Comments:
    #   Nodes and weights are originally taken from
    #   L. Zhang, T. Cui, and H. Liu, A set of symmetric quadrature rules
    #   on triangles and tetrahedra, J Comp Math (27) 2009
    #
    # References:
    #    Axel Grundmann, Michael Moeller,
    #    Invariant Integration Formulas for the N-Simplex by Combinatorial
    #    Methods, SIAM Journal on Numerical Analysis,
    #    Volume 15, Number 2, April 1978, pages 282-290.
    #    Arthur Stroud,
    #    Approximate Calculation of Multiple Integrals,
    #    Prentice Hall, 1971,
    #    Arthur H Stroud,
    #    A Fifth Degree Integration Formula for the N-Simplex,
    #    SIAM Journal on Numerical Analysis,
    #    Volume 6, Number 1, March 1969.
    #
    # Remark:
    #   This program is a supplement to the paper
    #   >> Efficient P1-FEM for any space dimension in Matlab <<
    #   by S. Beuter, and S. Funken. The reader should
    #   consult that paper for more information.
    #
    # Authors:
    #   S. Beuter, S. Funken 18-10-22

    maxOrder = 5
    if order > maxOrder
        println("Warning: Gauss type rule of order ",order," not supported. ")
        println("Using Gauss quadrature of order ",maxOrder," . ")
    end
    # Table of quadrature formulae on triangle of different order
    wdata = Array{Array}(undef,3)
    xdata = Array{Array}(undef,3)
    if order <= 1   # Reference: Grundmann&Moeller, Stroud (1971)
        wdata[1] = [1]
    elseif order <= 2   # Reference: Stroud (1971)
        xdata[2] = [1/(nD+1) * (1-1/sqrt(nD+2))]
        wdata[2] = [1/(nD+1)]
    elseif order <= 3     # Reference: Grundmann&Moeller, Stroud (1971)
        wdata[1] = [-0.25*(nD+1)*(nD+1)/(nD+2)]
        xdata[2] = [1/(nD+3)]
        wdata[2] = [0.25*(nD+3)*(nD+3)/((nD+1)*(nD+2))]
    elseif order <= 5   # Reference: Stroud (1969)
        wdata[1] = [(39*nD^2-307*nD+626)*(nD+1)^5/(72*prod(nD.+collect(1:5)))]
        xdata[2] = [((nD+4)-sqrt(15))/((nD+8)*nD+1);
            ((nD+4)+sqrt(15))/((nD+8)*nD+1)]
        wdata[2] = [(585-151*sqrt(15))*(93-7*nD+16*sqrt(15))*(nD+4+sqrt(15))^5;
                    (585+151*sqrt(15))*(93-7*nD-16*sqrt(15))*(nD+4-sqrt(15))^5] ./
                    (7560*prod(nD.+(1:5)))
        xdata[3] = [(nD+7+2*sqrt(15))/((nD+14)*nD-11);
            (nD+7-2*sqrt(15))/((nD+14)*nD-11)]
        wdata[3] = [(585+151*sqrt(15))*(nD+7-2*sqrt(15))^5;
            (585-151*sqrt(15))*(nD+7+2*sqrt(15))^5]./(1080*prod(nD.+(1:5)))
    end
    # Generate quadrature points and weights from data above
    lambda = zeros(0,nD+1)
    weights = []
    if isassigned(wdata,1)
        lambda = ones(1,nD+1)/(nD+1)
        weights = wdata[1]
    end
    if isassigned(wdata,2)
        perm = ones(Int64,nD+1).+Matrix(1I,nD+1,nD+1)
        points = [xdata[2] 1 .- nD .* xdata[2]]
        for i = 1:nD+1
            lambda = [lambda; points[:,perm[i,:]]]
            weights = [weights;wdata[2]]
        end
    end
    if isassigned(wdata,3)
        idx = reshape(hcat(vcat(1:nD.*(nD+1)÷2)*[1 1],
            vcat(collect(combinations(1:nD+1,2))'...)),:,2)
        perm = ones(Int64,maximum(idx[:,1]),maximum(idx[:,2]))
        for i = 1:size(idx,1)
            perm[idx[i,1],idx[i,2]] += 1
        end
        points = hcat(xdata[3],(1 .-(nD-1).*xdata[3])/2)
        for j = 1:size(perm,1)
            lambda = [lambda;points[:,perm[j,:]]]
            weights = [weights;wdata[3]]
        end
    end
    return lambda,weights
end
