include("infinite.jl")
include("../src/transfer_charge_labels.jl")
for path in ARGS
    started=time();actual=resolve_checkpoint(path);fingerprint=open(io->bytes2hex(sha256(io)),actual)
    psi=load_state(path);@assert hasqns(siteind(psi.AL,1)) "An original QN payload is required"
    meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    original=dag(input_inds(TransferMatrix(psi.AL)));n=length(psi.AL)
    A=InfiniteMPS([dense(psi.AL[j]) for j=1:n],translator(psi.AL));T=TransferMatrix(A)
    input=dag(input_inds(T));dims=Tuple(dim.(input));seed=7231
    start=ITensor(randn(MersenneTwister(seed),ComplexF64,dims...),input...)
    vals,vecs,info=KrylovKit.eigsolve(T,start,min(4,prod(dims)),:LM;tol=1e-9,krylovdim=30,maxiter=200)
    scale=val(qn(siteind(psi.AL,1),ITensors.Block(1)),"Sz")-val(qn(siteind(psi.AL,1),ITensors.Block(2)),"Sz")
    records=[]
    for k in sortperm(abs.(vals);rev=true)
        weights=transfer_charge_weights(vecs[k],original)
        push!(records,Dict("eigenvalue_real"=>real(vals[k]),"eigenvalue_imag"=>imag(vals[k]),"magnitude"=>abs(vals[k]),
            "residual"=>norm(T(vecs[k])-vals[k]*vecs[k]),"virtual_QN_flux_weights"=>Dict(string(q)=>p for (q,p) in weights),
            "virtual_number_difference_weights"=>Dict(string(q/scale)=>p for (q,p) in weights)))
    end
    @assert open(io->bytes2hex(sha256(io)),actual)==fingerprint
    output=Dict("source_checkpoint"=>actual,"source_payload_sha256"=>fingerprint,"family"=>meta["family"],"width"=>meta["width"],
        "cap"=>meta["cap"],"actual_chi"=>maxlinkdim(psi),"leading_modes"=>records,"converged_eigenpairs"=>info.converged,
        "physical_site_QN_up_minus_down"=>scale,"runtime_seconds"=>time()-started,
        "interpretation"=>"New eigenvector observable from unchanged AL tensors. Degenerate modes can mix U1 blocks; report their weights rather than assigning a unique charge. Virtual number differences do not classify microscopic CGS charges, operator overlaps, a physical gap, or the 2D phase.",
        "audit"=>run_provenance(seed=seed,solver="Library transfer eigensolve with original-QN eigenvector classification",
            settings=Dict("tol"=>1e-9,"krylovdim"=>30,"maxiter"=>200),initialization=path,conserved_quantum_numbers=["Original U1 virtual basis retained for classification"]))
    atomic_json(replace(path,".jls"=>"_transfer_charge_audit.json"),output)
    println("Classified ",length(records)," transfer modes from ",path);flush(stdout)
end
