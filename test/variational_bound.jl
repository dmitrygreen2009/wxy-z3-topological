using Test
isdefined(Main,:microscopic_mpo) || include("../src/model.jl")
@testset "Microscopic disjoint-dimer variational bound" begin
    # Each A star uses leg one and each B star leg two, all on distinct
    # physical gauges. Only those species-one exchanges have nonzero mean.
    s=siteinds("S=1/2",18;conserve_qns=true)
    pairs=[(1,13),(4,17),(7,16),(10,14)]
    pattern=fill("Dn",18)
    for (m,g) in pairs;pattern[m]="Up";end
    remaining=setdiff(1:18,vcat([collect(p) for p in pairs]...))
    for j in remaining[1:5];pattern[j]="Up";end
    psi=MPS(s,pattern)
    for (m,g) in pairs
        rotation=exp((pi/4)*(op("S-",s[m])*op("S+",s[g])-op("S+",s[m])*op("S-",s[g])))
        psi=apply(rotation,psi;cutoff=1e-14)
    end
    H=microscopic_mpo([[1,2,3],[1,5,3],[4,5,6],[4,2,6]],s)
    @test val(flux(psi),"Sz")==0
    @test norm(psi)≈1 atol=1e-12
    energy=real(inner(psi',H,psi))
    @test energy≈-4/sqrt(3) atol=1e-12
    atomic_json("results/variational_bound_audit.json",Dict("physical_spins"=>18,"nup"=>9,
        "bell_pairs_one_based_physical_indices"=>pairs,"endpoint_leg_table"=>[[1,2,3],[1,5,3],[4,5,6],[4,2,6]],
        "trial_energy"=>energy,"exact_trial_energy"=>-4/sqrt(3),"energy_error"=>abs(energy+4/sqrt(3)),
        "state_norm"=>norm(psi),"infinite_two_slice_bound"=>"E0_cell <= -4*width/sqrt(3)",
        "interpretation"=>"Half-filled product of distinct species-one Bell dimers; this is a variational upper bound, not an eigenstate or ground-state energy.",
        "audit"=>run_provenance(solver="ITensorMPS exact trial-state contractions",settings=Dict("gate_cutoff"=>1e-14,"energy_tolerance"=>1e-12),
            initialization="Deterministic QN product state and exact two-spin rotations",conserved_quantum_numbers=["N_up=9"])))
end
