# Serial recovery batch. Existing converged states and immutable checkpoints
# are retained; a failed sector does not prevent the other sectors being tried.
include("convert_winding_checkpoint.jl")
include("winding_qn_dmrg.jl")

function continue_winding_physical_sectors(;maxcap=512)
    outcomes=[]
    sources=Dict("zigzag"=>"results/zigzag_L2_w1_chi128_star.jls",
        "armchair"=>"results/armchair_L2_w1_chi256_star_Nup8.jls")
    for charge in 0:2,family in ["zigzag","armchair"]
        nup=family=="zigzag" ? 9 : 8
        path="results/$(family)_L2_w1_Nup$(nup)_winding$(charge)_qn_seed7254.json"
        status=Dict("family"=>family,"length"=>2,"width"=>1,"nup"=>nup,"winding_charge"=>charge)
        try
            prior=isfile(path) ? JSON3.read(read(path,String),Dict{String,Any}) : nothing
            if prior!==nothing
                @assert prior["family"]==family && prior["length"]==2 && prior["width"]==1 && prior["nup"]==nup
                @assert prior["audit"]["solver_settings"]["winding_charge"]==charge
            end
            if prior!==nothing && last(prior["records"])["converged_within_fixed_number_and_loop_sector"]
                # Reuse only a hash-checked physical-number/winding payload.
                load_valid_checkpoint(prior["checkpoint_file"])
                status["status"]="Previously converged state retained without optimization"
            else
                resume=prior===nothing ? nothing : prior["checkpoint_file"]
                if charge==0 && resume===nothing
                    resume=convert_winding_checkpoint(sources[family],0;target_nup=nup,maxdim=maxcap)
                elseif charge==0 && prior!==nothing && last(prior["records"])["energy_variance"]>1e-8
                    # A cold initialization can miss a CGS block. The validated
                    # physical ground checkpoint supplies a different exact seed.
                    resume=convert_winding_checkpoint(sources[family],0;target_nup=nup,maxdim=maxcap)
                end
                run_winding_qn(family,2,1,nup,charge;maxcap,resume,
                    noise=charge==0 ? [0.0] : [1e-5,1e-6,1e-7,0.0])
                status["status"]="Converged within fixed number and winding sector"
            end
        catch error
            status["status"]="Incomplete: checkpoints retained for targeted debugging"
            status["error"]=sprint(showerror,error,catch_backtrace())
            println("Sector incomplete: ",family," Nup=",nup," q=",charge," ",status["error"])
        end
        status["result_path"]=path;push!(outcomes,status)
        atomic_json("results/winding_physical_sector_recovery_batch.json",
            Dict("outcomes"=>outcomes,"git_commit"=>LAUNCH_REVISION,
                "interpretation"=>"Serial optimizer recovery audit; no MES, thermodynamic filling or phase certificate"))
        flush(stdout);GC.gc()
    end
    outcomes
end
if abspath(PROGRAM_FILE)==@__FILE__
    continue_winding_physical_sectors(;maxcap=isempty(ARGS) ? 512 : parse(Int,ARGS[1]))
end
