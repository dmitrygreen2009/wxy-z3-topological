# Complete finite virtual-space audit of the official transfer operator.
isdefined(@__MODULE__,:dense_transfer_matrix) || include("dense_transfer_matrix.jl")
function transfer_uniqueness_certificate(AL;maximum_dimension=256,tolerance=1e-10)
    T=TransferMatrix(AL);dimension=prod(dim.(dag(input_inds(T))))
    if dimension>maximum_dimension
        return Dict("status"=>"not computed: complete virtual spectrum exceeds resource limit",
            "transfer_dimension"=>dimension,"maximum_dimension"=>maximum_dimension,
            "fixed_point_uniqueness_certified"=>false,"primitive_peripheral_spectrum_certified"=>false)
    end
    M=dense_transfer_matrix(T;maximum_dimension)
    pairs=eigen(M);values=pairs.values
    residuals=[norm(M*pairs.vectors[:,j]-values[j]*pairs.vectors[:,j]) for j in eachindex(values)]
    k=argmin(abs.(values.-1));fixed_count=count(z->abs(z-1)<tolerance,values)
    peripheral_count=count(z->abs(abs(z)-1)<tolerance,values)
    reliable=maximum(residuals)<tolerance && abs(values[k]-1)<tolerance && maximum(abs.(values))<=1+tolerance
    Dict("status"=>"complete dense virtual spectrum","transfer_dimension"=>dimension,
        "maximum_dimension"=>maximum_dimension,"tolerance"=>tolerance,
        "values_real"=>real.(values),"values_imag"=>imag.(values),"eigenpair_residuals"=>residuals,
        "fixed_point_multiplicity_at_tolerance"=>fixed_count,"peripheral_multiplicity_at_tolerance"=>peripheral_count,
        "fixed_point_uniqueness_certified"=>reliable && fixed_count==1,
        "primitive_peripheral_spectrum_certified"=>reliable && fixed_count==1 && peripheral_count==1,
        "interpretation"=>"Complete finite virtual-space numerical certificate at the recorded tolerance; no physical bulk-gap or ground-state certificate")
end
