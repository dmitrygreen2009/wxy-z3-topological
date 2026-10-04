"""Library transfer spectrum with independent dominant-fixed-point starts.
A scalar Krylov space can hide repeated eigenvalues. Independent fixed points
protect xi extraction for noninjective states, such as a GHZ MPS.
"""
function transfer_spectrum(AL;tol=1e-9,seed=7103)
    T=TransferMatrix(AL)
    input=dag(input_inds(T));dimensions=Tuple(dim.(input))
    start(rngseed)=ITensor(randn(MersenneTwister(rngseed),ComplexF64,dimensions...),input...)
    count=min(4,prod(dimensions))
    values,vectors,info=KrylovKit.eigsolve(T,start(seed),count,:LM;tol,krylovdim=30,maxiter=200)
    permutation=sortperm(abs.(values);rev=true);values=values[permutation];vectors=vectors[permutation]
    residuals=[norm(T(vectors[k])-values[k]*vectors[k]) for k in eachindex(values)]
    fixed_vectors=[vectors[1]];fixed_values=[values[1]];fixed_residuals=[residuals[1]]
    if prod(dimensions)>1
        alt_values,alt_vectors,_=KrylovKit.eigsolve(T,start(seed+1),1,:LM;tol,krylovdim=30,maxiter=200)
        k=argmax(abs.(alt_values))
        alternative=alt_vectors[k];lambda=alt_values[k]
        alt_residual=norm(T(alternative)-lambda*alternative)
        overlap=abs(inner(vectors[1],alternative))/(norm(vectors[1])*norm(alternative))
        if abs(lambda-values[1])<10tol && overlap<1-1e-7 && max(residuals[1],alt_residual)<10tol
            push!(fixed_vectors,alternative);push!(fixed_values,lambda);push!(fixed_residuals,alt_residual)
            push!(values,lambda);push!(vectors,alternative);push!(residuals,alt_residual)
        end
    end
    gram=[inner(a,b)/(norm(a)*norm(b)) for a in fixed_vectors,b in fixed_vectors]
    permutation=sortperm(abs.(values);rev=true);values=values[permutation];residuals=residuals[permutation]
    ratio=length(values)>=2 ? abs(values[2]/values[1]) : 0.0
    xi=ratio==0 ? 0.0 : abs(ratio-1)<1e-12 ? Inf : -1/log(ratio)
    Dict("values"=>values,"residuals"=>residuals,"xi_cells"=>xi,"converged"=>info.converged,
        "fixed_point_rank_detected"=>length(fixed_vectors),"fixed_point_gram_eigenvalues"=>eigvals(Hermitian(gram)),
        "fixed_point_residuals"=>fixed_residuals,"subleading_ratio"=>ratio,
        "ratio_gap_above_residual_scale"=>length(values)>=2 && abs(1-ratio)>10maximum(residuals[1:2]))
end
