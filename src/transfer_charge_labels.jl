# Decompose a saved dense transfer eigenvector using the original QN basis.
# This is observable bookkeeping; all transfer contractions/eigensolves use
# ITensorInfiniteMPS/KrylovKit. It does not identify microscopic CGS charges.
function transfer_charge_weights(vector,qn_indices;name="Sz")
    @assert length(qn_indices)==2 && all(hasqns,qn_indices)
    labels=[]
    for i in qn_indices
        sign=ITensors.dir(i)==ITensors.Out ? 1 : -1
        push!(labels,vcat([fill(sign*val(qn(i,ITensors.Block(b)),name),ITensors.blockdim(i,ITensors.Block(b))) for b=1:ITensors.nblocks(i)]...))
    end
    dense_indices=ITensors.removeqns.(qn_indices)
    amplitudes=Array(vector,dense_indices...);total=sum(abs2,amplitudes)
    @assert total>0
    weights=Dict{Int,Float64}()
    for a in axes(amplitudes,1),b in axes(amplitudes,2)
        charge=labels[1][a]+labels[2][b]
        weights[charge]=get(weights,charge,0.)+abs2(amplitudes[a,b])/total
    end
    weights
end
