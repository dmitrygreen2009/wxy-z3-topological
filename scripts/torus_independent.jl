# Independent sparse ED. No ITensor, Python Hamiltonian, or benchmark vectors.
using LinearAlgebra, SparseArrays, Random, KrylovKit, JSON3
BLAS.set_num_threads(1);Random.seed!(20261003)
omega=exp(2pi*im/3)
coupling=[omega^((a-1)*(i-1))/sqrt(3) for a=1:3,i=1:3]
# Each physical edge is an endpoint pair, with a local leg at each end.
edges=[(1,1,2,1),(1,2,4,2),(1,3,2,3),
       (3,1,4,1),(3,2,2,2),(3,3,4,3)]
# Independently enumerate every full-Hilbert ket, using reversed bit order.
basis=[b for b=0:(2^18-1) if count_ones(b)==9]
@assert length(basis)==48620 && length(unique(basis))==48620
index=Dict(b=>j for (j,b) in enumerate(basis))
rows=Int[];cols=Int[];vals=ComplexF64[]
for (e,(v1,i1,v2,i2)) in enumerate(edges), (v,i) in [(v1,i1),(v2,i2)], a=1:3
    matter_site=3(v-1)+a;gauge_site=12+e
    m=18-matter_site;g=18-gauge_site
    for (j,b) in enumerate(basis)
        mu=(b>>m)&1;sigma=(b>>g)&1
        if mu==1 && sigma==0
            ket=b-(1<<m)+(1<<g)
            push!(rows,index[ket]);push!(cols,j);push!(vals,-coupling[a,i])
        elseif mu==0 && sigma==1
            ket=b+(1<<m)-(1<<g)
            push!(rows,index[ket]);push!(cols,j);push!(vals,-conj(coupling[a,i]))
        end
    end
end
H=sparse(rows,cols,vals,length(basis),length(basis))
@assert norm(H-H')<1e-12
start=KrylovKit.Block([randn(ComplexF64,length(basis)) for _=1:8])
energies,vectors,info=eigsolve(H,start,16,:SR;ishermitian=true,krylovdim=160,tol=1e-12,maxiter=200)
p=sortperm(energies);energies=energies[p];vectors=vectors[p]
V=hcat(vectors...)
residuals=[norm(H*vectors[k]-energies[k]*vectors[k]) for k in eachindex(energies)]
gram_error=norm(V'*V-I)
@assert info.converged>=16
@assert maximum(residuals[1:16])<1e-9 && gram_error<1e-8
reference=vcat(fill(-7.364534839426,4),[-7.342883075944,-7.337288472473],fill(-7.296046028855,4))
@assert maximum(abs.(energies[1:10]-reference))<1e-9
r=Dict("solver"=>"Independent Julia/KrylovKit eight-vector BlockLanczos", "dimension"=>length(basis),
    "energies"=>energies,"residuals"=>residuals,"gram_error"=>gram_error,"converged"=>info.converged,
    "corrected_benchmark_first10"=>reference,"corrected_benchmark_error"=>maximum(abs.(energies[1:10]-reference)))
open("results/torus_independent.json","w") do io;JSON3.write(io,r);end
println(JSON3.write(r));flush(stdout)
