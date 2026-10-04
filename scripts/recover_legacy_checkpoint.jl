isdefined(Main,:cylinder) || include("../src/model.jl")
using SHA
function recover_zigzag_reference()
    legacy="results/zigzag_L4_w1_chi256_star.jls"
    try
        load_state(legacy)
        return
    catch error
        reason=sprint(showerror,error)
        replacement="results/zigzag_L4_w1_chi256_star_seed7104.jls"
        good=JSON3.read(read(replace(replacement,".jls"=>".json"),String),Dict{String,Any})
        psi=load_state(replacement);lat=cylinder("zigzag",4,1;ordering="star")
        @assert length(psi)==lat.n && val(flux(psi),"Sz")==2cld(lat.n,2)-lat.n
        H=MPO(lat.os,siteinds(psi));E=real(inner(psi',H,psi));S,p=entropy_at(psi,lat.cut)
        @assert abs(E-good["records"][end]["energy"])<1e-10
        @assert abs(S-good["records"][end]["entropy"])<1e-10
        archive="results/legacy_invalid_checkpoint_audit";mkpath(archive)
        original_json=joinpath(archive,"zigzag_L4_w1_chi256_star_original_result.json")
        isfile(original_json) || cp(replace(legacy,".jls"=>".json"),original_json)
        original_binary="results/checkpoints/legacy_invalid_zigzag_L4_w1_chi256_star.jls"
        bytes=isfile(legacy) ? filesize(legacy) : nothing
        digest=isfile(legacy) ? open(io->bytes2hex(sha256(io)),legacy) : nothing
        isfile(legacy) && mv(legacy,original_binary;force=false)
        recovery=Dict("error"=>reason,"unreadable_original_bytes"=>bytes,"unreadable_original_sha256"=>digest,
            "original_payload_preserved_at"=>original_binary,"original_result_preserved_at"=>original_json,
            "replacement_source"=>replacement,"replacement_seed"=>7104,"replacement_energy"=>E,"replacement_entropy"=>S,
            "git_commit"=>LAUNCH_REVISION,"interpretation"=>"Unreadable under the pinned working environment; cause not inferred. Original bytes and scalar results preserved. Alias now resolves a independently validated seed7104 state, not the original seed7103 wavefunction.")
        good["legacy_checkpoint_recovery"]=recovery
        atomic_json(replace(legacy,".jls"=>".json"),good)
        atomic_json(joinpath(archive,"recovery.json"),recovery)
        println("Recovered legacy alias from validated seed7104 state; exact prior error: ",reason);flush(stdout)
    end
end
