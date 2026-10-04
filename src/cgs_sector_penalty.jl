# Exact commuting projector penalty. All MPO operations use ITensorMPS.
isdefined(@__MODULE__, :cgs_product_mpo) || include("cgs_operator_mpo.jl")
function cgs_full_spec(lat,p)
    @assert length(p)==lat.n-3lat.nv
    @assert all(mod(sum(p[e] for e in es),3)==0 for es in lat.legs)
    ip=invperm(lat.order)
    triplets=Tuple{Int,Int,Int}[];gauges=Tuple{Int,Int}[]
    for v=1:lat.nv
        c=mod(p[lat.legs[v][1]],3);q=mod(p[lat.legs[v][2]]-c,3)
        c==0 && q==0 && continue
        js=[ip[3(v-1)+a] for a=1:3]
        @assert js==collect(first(js):first(js)+2) "Matter blocks must be contiguous and species ordered"
        push!(triplets,(first(js),c,q))
    end
    for e in eachindex(p)
        mod(p[e],3)==0 || push!(gauges,(ip[3lat.nv+e],mod(p[e],3)))
    end
    (;triplets,gauges,first_site=1,last_site=lat.n)
end
function cgs_sector_penalty(lat,sites,p;charge,strength=nothing,cutoff=1e-14)
    @assert charge in 0:2
    # Each Hermitian exchange has norm |W_ai|. B bounds ||H||.
    bound=lat.nv*sum(abs,W)
    lambda=strength===nothing ? 2bound+1 : Float64(strength)
    @assert lambda>2bound "Penalty must separate all sectors rigorously"
    U=cgs_product_mpo(sites,cgs_full_spec(lat,p);cutoff)
    U2=cgs_product_mpo(sites,cgs_full_spec(lat,mod.(2p,3));cutoff)
    identity=cgs_product_mpo(sites,(;triplets=Tuple{Int,Int,Int}[],gauges=Tuple{Int,Int}[],first_site=1,last_site=lat.n);cutoff)
    z=cis(2pi*charge/3)
    a=deepcopy(U);b=deepcopy(U2)
    a[1]*=-lambda*conj(z)/3;b[1]*=-lambda*conj(z)^2/3
    identity[1]*=2lambda/3
    (;terms=[identity,a,b],U,U2,strength=lambda,hamiltonian_norm_bound=bound,charge,
        interpretation="H + lambda(I-P_q); exactly H within sector q. MPO updates need not stay pure; final leakage and stationarity must be validated.")
end
