include("infinite.jl")
function run_infinite_qn(family,w,chi)
    Random.seed!(7103)
    n=18w
    function initstate(j)
        k=mod(j-1,9)+1
        k<=3 ? "Up" : 6<=k<=8 ? "Dn" : k==4 ? "Dn" : k==5 ? "Up" : (iseven(fld(j-1,9w)) ? "Dn" : "Up")
    end
    @assert count(j->initstate(j)=="Up",1:n)==n÷2
    s=infsiteinds("S=1/2",n;conserve_qns=true,initstate)
    psi=InfMPS(s,initstate)
    H=InfiniteSum{MPO}(infinite_opsum(family,w;ordering="star"),s)
    for (iteration,cap) in enumerate(unique(vcat([min(c,chi) for c in [2,4,8,16,32,64,128]],chi)))
        estimate_memory(psi,cap;label="$(family) infinite w$(w) expansion")
        psi=subspace_expansion(psi,H;cutoff=1e-10,maxdim=cap)
        psi=audited_vumps(H,psi;family,w,cap,tag="_qn_star_active",ordering="star",tol=1e-7,maxiter=30,
            solver_tol=x->1e-10)
        measure_infinite(psi,H,family,w,cap,iteration;tag="_qn_star_active",ordering="star")
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;run_infinite_qn(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]));end
