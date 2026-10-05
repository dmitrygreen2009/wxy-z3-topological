# Small algebra/update fixture, not a converged cylinder calculation.
using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
isdefined(Main,:audited_vumps) || include("../scripts/infinite.jl")
isdefined(Main,:infinite_cell_product) || include("../src/infinite_block_initializer.jl")
include("../src/infinite_winding_basis.jl")
fixture_records=[]

@testset "Individual infinite microscopic winding-cell constraint" begin
    for family in ["zigzag","armchair"],q in 0:1
        charges=family=="zigzag" ? [q,q] : [q]
        phi,layout=infinite_winding_initial_cell(family,1,6,charges;seed=7259,linkdims=2)
        @test iszero(flux(phi))
        @test (18+sum(expect(phi,"Sz"))*2)/2≈6 atol=1e-10
        sites=siteinds(phi);Hrot=MPO(infinite_real_charge_opsum(family,1;periodic=true),sites)
        physical=dense(phi);ss=siteinds(physical);gates=ITensor[]
        B=matter_charge_basis()[8:-1:1,8:-1:1]
        for row=0:1,first in [9row+1,9row+6]
            local_sites=ss[first:first+2]
            push!(gates,ITensor(reshape(B,2,2,2,2,2,2),prime.(local_sites)...,dag.(local_sites)...))
        end
        physical=apply(gates,physical;cutoff=1e-14,maxdim=256)
        Hold=MPO(infinite_opsum(family,1;periodic=true,ordering="star"),ss)
        @test inner(phi',Hrot,phi)≈inner(physical',Hold,physical) atol=1e-10
        rotated_density=expect(phi,"Sz");physical_density=expect(physical,"Sz")
        for row=0:1,first in [9row+1,9row+6]
            for a=0:2
                @test physical_density[first+a]≈sum(rotated_density[first:first+2])/3 atol=1e-10
            end
        end
        for group=1:length(charges)
            offset=family=="zigzag" ? (group-1)*9 : 0
            spec=merge(cgs_circumference_spec(family,1;offset),(;first_site=1,last_site=18))
            originalU=cgs_product_mpo(ss,spec)
            @test inner(physical',originalU,physical)≈cis(2pi*q/3) atol=1e-10
        end
        psi=infinite_cell_product(phi)
        @test !isempty(infinite_winding_boundary_support(psi,layout))
        @test infinite_number_background(psi)["physical_number_density"]≈1/3
        function check_loops(state,stage)
            s=siteinds(only,state)
            measurements=[]
            for group=1:length(charges),offset in [0,layout.n]
                z=expect(state,cgs_product_mpo(s,infinite_winding_diagonal_spec(layout,group;offset)))
                @test z≈cis(2pi*q/3) atol=1e-10
                @test abs(1-abs2(z))<1e-10
                push!(measurements,Dict("group"=>group,"translated_sites"=>offset,"real"=>real(z),"imag"=>imag(z),"variance"=>1-abs2(z)))
            end
            push!(fixture_records,Dict("family"=>family,"width"=>1,"cell_spins"=>18,"cell_slices"=>2,
                "target_loop_charges"=>charges,"stage"=>stage,"measurements"=>measurements,
                "boundary_support"=>infinite_winding_boundary_support(state,layout),
                "number_background"=>infinite_number_background(state),"physical_nup_per_cell"=>6,
                "physical_sz_per_cell"=>-3,"filling_fraction"=>1/3,"filling_selection"=>"Explicit fixed U1 background fixture, not ground selection"))
        end
        check_loops(psi,"canonical product of exact finite cells")
        H=InfiniteSum{MPO}(infinite_real_charge_opsum(family,1),siteinds(only,psi))
        psi=subspace_expansion(psi,H;cutoff=1e-12,maxdim=3)
        @test !isempty(infinite_winding_boundary_support(psi,layout))
        check_loops(psi,"official subspace expansion, maxdim=3")
        # One package update checks symmetry closure, not variational convergence.
        epsL=fill(1e-7,layout.n);epsR=fill(1e-7,layout.n)
        psi,_=ITensorInfiniteMPS.tdvp_iteration(ITensorInfiniteMPS.vumps_solver,H,psi;
            (ϵᴸ!)=epsL,(ϵᴿ!)=epsR,multisite_update_alg="sequential",solver_tol=x->1e-10,time_step=-Inf,eager=true)
        @test !isempty(infinite_winding_boundary_support(psi,layout))
        check_loops(psi,"one official sequential infinite update")
        @test infinite_number_background(psi)["physical_number_density"]≈1/3
    end
end
fixture_output=Dict("records"=>fixture_records,
    "audit"=>run_provenance(;seed=7259,solver="Official infinite canonicalization/expansion/one-update symmetry fixture",
        settings=Dict("cutoff"=>1e-12,"maxdim"=>3,"iterations"=>1,"local_tolerance"=>1e-10,"loop_tolerance"=>1e-10),
        initialization="Random finite cells with exact physical-number/winding QNs",conserved_quantum_numbers=["Physical mean N_up=6 per 18-site cell","Individual disjoint microscopic winding charges"]),
    "interpretation"=>"Algebra and package-update closure only; not a converged cylinder, selected thermodynamic filling, MES or phase certificate.")
atomic_json("results/infinite_winding_basis_fixture.json",fixture_output)
println("INFINITE_WINDING_FIXTURE_JSON=",JSON3.write(fixture_output))
