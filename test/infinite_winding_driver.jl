# One tiny integration fixture checks hooks, basis interpretation and resume data.
# It deliberately cannot qualify as a converged physical cylinder result.
using Test
include("../scripts/infinite_winding_qn.jl")
@testset "Constrained infinite driver audit integration fixture" begin
    records=run_infinite_winding_qn("zigzag",1,6,0,2;maxiter=1,validation_fixture=true)
    @test length(records)==1
    @test records[1]["cap"]==2
    checkpoint=records[1]["resume_checkpoint"]
    meta=JSON3.read(read(replace(checkpoint,".jls"=>".json"),String),Dict{String,Any})
    @test checkpoint_basis(meta)=="exact_matter_charge_basis"
    @test meta["number_background"]["target_average_nup_per_cell"]==6
    @test haskey(meta["iterations"][1],"additional_symmetry_validation")
    @test meta["solver_settings"]["individual_loop_charges"]==[0,0]
    @test hasqns(siteind(load_state(checkpoint).AL,1))
end
