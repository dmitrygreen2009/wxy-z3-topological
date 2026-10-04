include("cylinders.jl")
include("continue.jl")
reference=JSON3.read(read("results/small_cylinders_ed.json",String))
for family in ["zigzag","armchair"]
    stem="results/$(family)_L2_w1_chi128_star"
    isfile(stem*".jls") || run_cylinder(family,2,1,128;ordering="star")
    for cap in [128,256,512]
        if cap>128
            continue_cylinder(stem*".jls",cap)
            stem="results/$(family)_L2_w1_chi$(cap)_star"
        end
        meta=JSON3.read(read(stem*".json",String),Dict{String,Any})
        psi=load_state(stem*".jls")
        H=MPO(cylinder(family,2,1;ordering="star").os,siteinds(psi))
        energy=real(inner(psi',H,psi))
        variance=real(inner(H,psi,H,psi))-energy^2
        error=abs(energy-reference[family]["energies"][1])
        meta["independent_ed_error"]=error;meta["energy_variance"]=variance
        meta["independent_ed_passed"]=error<1e-7 && abs(variance)<1e-6
        open(stem*".json","w") do io;JSON3.write(io,meta);end
        if meta["independent_ed_passed"]
            println("VALIDATED SMALL CYLINDER ",family," chi=",cap," error=",error," variance=",variance)
            flush(stdout);break
        end
        @assert cap<512 "Small-cylinder ED comparison still fails at chi=512"
    end
end
