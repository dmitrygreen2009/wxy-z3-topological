include("infinite.jl")
include("../src/dense_transfer_matrix.jl")
for path in ARGS
 started=time();actual=resolve_checkpoint(path);source_sha=open(io->bytes2hex(sha256(io)),actual)
 psi=load_state(path);n=length(psi.AL)
 A=InfiniteMPS([dense(psi.AL[j]) for j=1:n],translator(psi.AL));T=TransferMatrix(A)
 matrix=dense_transfer_matrix(T);dimension=size(matrix,1)
 eigenpairs=eigen(matrix);permutation=sortperm(abs.(eigenpairs.values);rev=true)
 values=eigenpairs.values[permutation];vectors=eigenpairs.vectors[:,permutation]
 residuals=[norm(matrix*vectors[:,k]-values[k]*vectors[:,k])/norm(vectors[:,k]) for k=1:dimension]
 ratio=dimension>1 ? abs(values[2]/values[1]) : 0.
 xi=ratio==0 ? 0. : abs(ratio-1)<1e-12 ? Inf : -1/log(ratio)
 meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
 slices=n÷(9meta["width"])
 r=Dict("source_checkpoint"=>actual,"source_payload_sha256"=>source_sha,"family"=>meta["family"],"width"=>meta["width"],
 "cell_slices"=>slices,"transfer_dimension"=>dimension,"eigenvalues_real"=>real.(values),"eigenvalues_imag"=>imag.(values),
 "eigenvalue_residuals"=>residuals,"maximum_eigenpair_residual"=>maximum(residuals),
 "full_peripheral_magnitude_count_at_1e-9"=>count(z->abs(abs(z/values[1])-1)<1e-9,values),
 "xi_cells"=>isfinite(xi) ? xi : "Inf","xi_slices"=>isfinite(xi) ? slices*xi : "Inf","runtime_seconds"=>time()-started,
 "interpretation"=>"Full finite virtual-space spectrum of the unchanged AL state via official TransferMatrix actions and LinearAlgebra.eigen. Numerical peripheral count at the stated tolerance; no physical bulk gap or ground-state certificate.",
 "audit"=>run_provenance(solver="Independent dense virtual transfer eigenspectrum",initialization=actual,settings=Dict("maximum_virtual_dimension"=>1024,"peripheral_magnitude_tolerance"=>1e-9),conserved_quantum_numbers=[]))
 @assert open(io->bytes2hex(sha256(io)),actual)==source_sha
 atomic_json(replace(path,".jls"=>"_dense_transfer_audit.json"),r)
 println("DENSE TRANSFER dimension=",dimension," peripheral_count=",r["full_peripheral_magnitude_count_at_1e-9"]," xi_cells=",r["xi_cells"]);flush(stdout)
end
