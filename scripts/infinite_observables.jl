include("infinite.jl")
include("../src/infinite_observables.jl")
function measure_infinite_observables(path)
    started=time();meta=JSON3.read(read(replace(path,".jls"=>".json"),String),Dict{String,Any})
    @assert get(meta,"infinite_ordering","matter_first")=="star" "Observable locations require the recorded star ordering"
    actual=resolve_checkpoint(path);source_sha=open(io->bytes2hex(sha256(io)),actual)
    psi=load_state(actual);psi=InfiniteCanonicalMPS(dense(psi.AL),dense(psi.C),dense(psi.AR))
    s=siteinds(only,psi);n=length(s);w=meta["width"]
    loops=meta["family"]=="zigzag" && w==1 ? [expect(psi,infinite_twocycle_mpo(s,base)) for base=0:9:n-1] : ComplexF64[]
    result=Dict("family"=>meta["family"],"width"=>w,"cell_slices"=>n÷(9w),"chi"=>meta["chi"],"source_file"=>path,
        "source_payload"=>actual,"source_payload_sha256"=>source_sha,
        "parallel_two_cycle_expectations_real"=>real.(loops),"parallel_two_cycle_expectations_imag"=>imag.(loops),
        "original_u1_conserving_ansatz"=>get(meta,"u1_conserving_ansatz",nothing),"cycle_charge_probabilities"=>[[real((1+2real(conj(cis(2pi*q/3))*z))/3) for q=0:2] for z in loops],"source_canonical_error"=>meta["canonical_error"],"source_solver_residual"=>get(meta,"solver_residual",nothing),
        "runtime_seconds"=>time()-started,"audit"=>run_provenance(;seed=nothing,solver="ITensorInfiniteMPS library observable contractions",settings=Dict("MPO_cutoff"=>1e-14),initialization="Measurement-only dense copy of $path",conserved_quantum_numbers=[]),
        "interpretation"=>"Microscopic local loop charges of the saved variational state. Zigzag width one has parallel-edge local S3 symmetry. Unconverged variational states do not establish a thermodynamic phase.")
    @assert open(io->bytes2hex(sha256(io)),actual)==source_sha "Source checkpoint changed during measurement"
    atomic_json(replace(path,".jls"=>"_observables.json"),result)
    println(meta["family"]," w=",w," chi=",meta["chi"]," loops=",loops);flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__;for path in ARGS;measure_infinite_observables(path);end;end
