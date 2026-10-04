# Independent finite-dimensional audit of the official transfer operator.
# This forms a matrix by applying the existing ITensorInfiniteMPS operator;
# it does not implement tensor-network contraction or optimization.
function dense_transfer_matrix(T;maximum_dimension=1024)
    ii=dag(input_inds(T));dims=Tuple(dim.(ii));dimension=prod(dims)
    @assert dimension<=maximum_dimension "Dense transfer audit exceeds configured virtual-space limit"
    matrix=zeros(ComplexF64,dimension,dimension)
    for k=1:dimension
        data=zeros(ComplexF64,dimension);data[k]=1
        basis=ITensor(reshape(data,dims),ii...)
        matrix[:,k]=vec(Array(T(basis),ii...))
    end
    return matrix
end
