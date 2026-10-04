include("cylinders.jl")
include("continue.jl")
function resource_checked(label,calculation)
    try
        calculation()
    catch error
        if error isa ResourceLimitError
            record=Dict("label"=>label,"status"=>"resource_limited","error"=>sprint(showerror,error),
                "git_commit"=>LAUNCH_REVISION,"recorded_utc"=>string(now(UTC)),"existing_results_preserved"=>true)
            open("results/resource_limited_points.jsonl","a") do io;println(io,JSON3.write(record));end
            println(JSON3.write(record));flush(stdout)
        else
            rethrow()
        end
    end
end
# Common chi=512 for all four circumferences; use completed lower-chi states.
for family in ["zigzag","armchair"]
    resource_checked("$(family) L4 w4 chi512",()->continue_cylinder("results/$(family)_L4_w4_chi256_star.jls",512))
end
# Longitudinal checks where widths three and four are still short cylinders.
for family in ["zigzag","armchair"],(L,w) in [(8,3),(6,4),(8,4)]
    resource_checked("$(family) L$(L) w$(w) chi256",()->run_cylinder(family,L,w,256;ordering="star"))
    point="results/$(family)_L$(L)_w$(w)_chi256_star.json"
    isfile(point) || continue
    resource_checked("$(family) L$(L) w$(w) chi512",()->continue_cylinder(replace(point,".json"=>".jls"),512))
end
