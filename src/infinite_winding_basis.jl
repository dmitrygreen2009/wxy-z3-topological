# Exact microscopic representation and initialization; official libraries
# construct/contract every MPS/MPO and perform all infinite updates.
isdefined(@__MODULE__,:cgs_circumference_spec) || include("cgs_operator_mpo.jl")
isdefined(@__MODULE__,:winding_initial_state) || include("matter_charge_basis.jl")

function infinite_real_charge_opsum(family,w;cell_slices=2,periodic=false)
    n=9w*cell_slices;os=OpSum()
    a=(1-inv(sqrt(3)))/2;b=2*(1+inv(sqrt(3)));d=2/sqrt(3)
    index(j)=periodic ? mod1(j,n) : j
    for p in infinite_pairs(family,w;ordering="star",cell_slices)
        p.a==1 || continue
        i=p.i;other=filter(!=(i),collect(1:3))
        si=index(p.m+i-1);sj=index(p.m+other[1]-1);sk=index(p.m+other[2]-1);sg=index(p.g)
        add!(os,-a,"S-",si,"S+",sg);add!(os,-a,"S+",si,"S-",sg)
        add!(os,-b,"S-",si,"Sz",sj,"Sz",sk,"S+",sg)
        add!(os,-b,"S+",si,"Sz",sj,"Sz",sk,"S-",sg)
        add!(os,-d,"S+",si,"S-",sj,"S-",sk,"S+",sg)
        add!(os,-d,"S-",si,"S+",sj,"S+",sk,"S-",sg)
    end
    if !periodic
        for j=1:n;add!(os,0.,"Id",j);end
    end
    os
end

function infinite_winding_cell_layout(family,w,cell_slices,charges)
    step=family=="zigzag" ? 1 : 2
    @assert family in ["zigzag","armchair"] && cell_slices%step==0
    groups=cell_slices÷step
    @assert length(charges)==groups && all(q in 0:2 for q in charges)
    @assert groups<=3 "Pinned ITensor QNs allow at most four fields including physical number; enlarge support through a separately validated alternative if needed"
    n=9w*cell_slices;weights=zeros(Int,n,groups);background=zeros(Int,n,groups)
    ends=Int[]
    for g=1:groups
        offset=(g-1)*step*9w;spec=cgs_circumference_spec(family,w;offset)
        for (first,c,q) in spec.triplets,a=0:2;weights[first+a,g]=mod(c+q*a,3);end
        for (j,p) in spec.gauges;weights[j,g]=mod(p,3);end
        last=g*step*9w;push!(ends,last);background[last,g]=charges[g]
        @assert weights[last,g]==0
    end
    (;weights,background,ends,charges=collect(charges),field_names=["Wind$g" for g=1:groups],n)
end

function infinite_winding_initial_cell(family,w,nup,charges;cell_slices=2,seed=7259,linkdims=2)
    Random.seed!(seed);layout=infinite_winding_cell_layout(family,w,cell_slices,charges)
    @assert 0<nup<layout.n
    # Initial allocation is not additional conserved per-loop physical number.
    groups=length(charges);width=layout.n÷groups
    counts=fill(nup÷groups,groups)
    for g=1:nup%groups;counts[g]+=1;end
    state=String[]
    for g=1:groups
        append!(state,winding_initial_state(layout.weights[(g-1)*width+1:g*width,g],counts[g],charges[g]))
    end
    original=siteinds("S=1/2",layout.n;conserve_qns=true)
    shifted=ITensorInfiniteMPS.shift_flux_to_zero(original,j->state[j])
    sites=Index[]
    for j=1:layout.n
        blocks=Pair{QN,Int}[]
        for bit=1:-1:0
            sz=val(qn(shifted[j],ITensors.Block(bit==1 ? 1 : 2)),"Sz")
            fields=Tuple[("Sz",sz)]
            for g=1:groups
                push!(fields,(layout.field_names[g],mod(layout.weights[j,g]*bit-layout.background[j,g],3),3))
            end
            push!(blocks,QN(fields...)=>1)
        end
        push!(sites,Index(blocks,tags(shifted[j])))
    end
    sites=collect(ITensorInfiniteMPS.infsiteinds(sites))
    phi=random_mps(Float64,sites,state;linkdims)
    @assert iszero(flux(phi))
    phi,layout
end

function infinite_winding_boundary_support(psi,layout)
    records=[]
    for tensors in [psi.AL,psi.AR],j in layout.ends
        link=commonind(tensors[j],tensors[j+1]);@assert hasqns(link)
        charges=Dict(name=>sort(unique(mod(val(qn(link,ITensors.Block(b)),name),3) for b=1:nblocks(link))) for name in layout.field_names)
        @assert all(support==[0] for support in Base.values(charges)) "Individual winding-boundary support expanded; do not claim sector preservation"
        push!(records,Dict("bond"=>j,"dimension"=>dim(link),"charge_support"=>charges))
    end
    records
end

function infinite_winding_diagonal_spec(layout,group;offset=0)
    active=findall(!iszero,layout.weights[:,group])
    (;triplets=Tuple{Int,Int,Int}[],gauges=[(j+offset,layout.weights[j,group]) for j in active],
        first_site=first(active)+offset,last_site=last(active)+offset)
end
