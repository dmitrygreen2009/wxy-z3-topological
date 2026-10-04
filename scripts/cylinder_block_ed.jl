# Independent full-spin ED from the machine-readable endpoint records.
# Block starts retain degeneracies that scalar Lanczos may omit.
using LinearAlgebra, SparseArrays, Random, KrylovKit, JSON3, SHA
BLAS.set_num_threads(1)
function cylinder_block_ed(family;nup_override=nothing)
    started=time();seed=7115;Random.seed!(seed)
    path="geometry/$(family)_L2_w1_star.json"
    geometry=JSON3.read(read(path,String),Dict{String,Any})
    n=geometry["physical_spins"];nup=nup_override===nothing ? cld(n,2) : nup_override
    @assert n<=22
    basis=[b for b=0:(2^n-1) if count_ones(b)==nup]
    lookup=Dict(b=>k for (k,b) in enumerate(basis))
    matter=Dict((m["vertex_id"],m["species_a"])=>m["physical_spin_id"] for m in geometry["matter_spins"])
    rows=Int[];cols=Int[];data=ComplexF64[]
    omega=exp(2pi*im/3)
    for edge in geometry["shared_gauge_edges"], endpoint in edge["endpoints"], a=1:3
        v=endpoint["vertex_id"];i=endpoint["local_leg_i"]
        @assert endpoint["W_convention"]=="W[a,i]"
        # Reverse bit order independently of the Python Hamiltonian.
        m=n-matter[(v,a)];g=n-edge["physical_spin_id"]
        coupling=omega^((a-1)*(i-1))/sqrt(3)
        for (column,b) in enumerate(basis)
            mu=(b>>m)&1;sigma=(b>>g)&1
            mu==sigma && continue
            push!(rows,lookup[xor(b,(1<<m)|(1<<g))]);push!(cols,column)
            push!(data,mu==1 ? -coupling : -conj(coupling))
        end
    end
    H=sparse(rows,cols,data,length(basis),length(basis))
    herm=norm(H-H');@assert herm<1e-12
    empty!(rows);empty!(cols);empty!(data);GC.gc()
    initial=KrylovKit.Block([randn(ComplexF64,length(basis)) for _=1:8])
    values,vectors,info=eigsolve(H,initial,16,:SR;ishermitian=true,krylovdim=160,tol=1e-12,maxiter=200)
    permutation=sortperm(values);values=values[permutation];vectors=vectors[permutation]
    residuals=[norm(H*v-values[k]*v) for (k,v) in enumerate(vectors)]
    gram=norm(hcat(vectors...)'*hcat(vectors...)-I)
    @assert info.converged>=16 && maximum(residuals)<1e-9 && gram<1e-8
    reference=if nup_override===nothing
        JSON3.read(read("results/small_cylinders_ed.json",String))[family]["energies"][1]
    else
        scan=JSON3.read(read("results/small_cylinder_number_sectors.json",String))[family*"_L2_w1"]["points"]
        only(p["energy"] for p in scan if p["nup"]==nup)
    end
    @assert abs(values[1]-reference)<1e-9
    ground_count=count(e->abs(e-values[1])<1e-9,values)
    @assert ground_count<16
    result=Dict("family"=>family,"length"=>2,"width"=>1,"physical_spins"=>n,"nup"=>nup,
        "dimension"=>length(basis),"energies"=>values,"residuals"=>residuals,"gram_error"=>gram,
        "ground_multiplicity"=>ground_count,"gap_within_number_sector"=>values[ground_count+1]-values[1],
        "hermiticity_error"=>herm,"runtime_seconds"=>time()-started,
        "solver"=>"Independent Julia/KrylovKit eight-vector BlockLanczos","random_seed"=>seed,
        "solver_settings"=>Dict("requested_eigenpairs"=>16,"krylovdim"=>160,"tol"=>1e-12,"maxiter"=>200),
        "initialization"=>"Eight independent complex Gaussian block vectors","conserved_quantum_numbers"=>["N_up=$nup"],
        "git_commit"=>strip(read(`git rev-parse HEAD`,String)),"julia_version"=>string(VERSION),
        "geometry_file"=>path,"geometry_sha256"=>bytes2hex(sha256(read(path))),
        "source_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "interpretation"=>"Finite open-cylinder spectrum only; last returned excited cluster may be incomplete. No thermodynamic phase inference.")
    sector_suffix=nup_override===nothing ? "" : "_Nup$(nup)"
    open("results/$(family)_L2_w1$(sector_suffix)_block_ed.json","w") do io;JSON3.write(io,result);end
    println(family," ground multiplicity=",ground_count," gap=",result["gap_within_number_sector"]," max residual=",maximum(residuals));flush(stdout)
end
for argument in ARGS
    parts=split(argument,":");cylinder_block_ed(parts[1];nup_override=length(parts)>1 ? parse(Int,parts[2]) : nothing);GC.gc()
end
