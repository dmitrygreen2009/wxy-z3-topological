include("../src/model.jl")
Random.seed!(7103)
BLAS.set_num_threads(1)
cases=[("star",[[1,2,3]],3,-2.403092127540),
       ("pair",[[1,2,3],[1,4,5]],5,-4.3967655435),
       ("torus",[[1,2,3],[1,5,3],[4,5,6],[4,2,6]],9,-7.364534839425767)]
out=Dict()
for (name,legs,nup,target) in cases
    n=3length(legs)+maximum(vcat(legs...))
    sites=siteinds("S=1/2",n;conserve_qns=true)
    H=microscopic_mpo(legs,sites)
    state=[i<=nup ? "Up" : "Dn" for i=1:n]
    psi=random_mps(sites,state;linkdims=8)
    e,psi=dmrg(H,psi;nsweeps=24,maxdim=[16,32,64,128,256,512],cutoff=1e-13,
        noise=[1e-5,1e-6,1e-7,0],eigsolve_krylovdim=12,outputlevel=1)
    variance=real(inner(H,psi,H,psi)-e^2)
    out[name]=Dict("energy"=>e,"target_error"=>abs(e-target),"variance"=>variance,"bond_dimension"=>maxlinkdim(psi))
    out[name]["verification_scope"]="ITensor ground energy and variance only; the full torus spectrum hard gate is scripts/torus_gate.py"
    open("results/itensor_validation.json","w") do io; JSON3.write(io,out); end
    @assert abs(e-target)<1e-9 "$name ground energy mismatch"
    @assert abs(variance)<1e-8
end
