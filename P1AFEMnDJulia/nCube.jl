function nCube(n,s,pattern)
    coordinates, elements = nSubCube(n,s,pattern)
    boundary = zeros(0,n)
    for k = 1:n+1
        bdry = elements[:,[1:k-1;k+1:n+1]]
        midpts = zeros(size(bdry,1),n)
        for j=1:n
            midpts = midpts + coordinates[bdry[:,j],:]
        end
        midpts = midpts./n
        idx = (maximum(abs.(midpts),dims=2) .>= 1-1e-4) .| (minimum(midpts,dims=2) .<= 1e-4)
        boundary = [boundary; bdry[idx[:],:]]
    end
    return (coordinates,elements,boundary)
end
