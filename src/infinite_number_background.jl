# The library may shift and scale physical site QNs to represent a rational
# filling. Zero virtual flux therefore does not by itself mean physical Sz=0.
function infinite_number_background(psi)
    s=siteinds(only,psi);n=length(s)
    if !hasqns(s[1])
        return Dict("label"=>"unrestricted","density_enforced"=>false)
    end
    up=val(qn(s[1],ITensors.Block(1)),"Sz");down=val(qn(s[1],ITensors.Block(2)),"Sz")
    @assert up>down
    @assert all(val(qn(j,ITensors.Block(1)),"Sz")==up && val(qn(j,ITensors.Block(2)),"Sz")==down for j in s[Cell(1)])
    density=(-down)//(up-down);number=n*density
    @assert denominator(number)==1
    label=2numerator(number)==n ? "Sz0" : "rho$(numerator(density))of$(denominator(density))"
    Dict("label"=>label,"density_enforced"=>true,"physical_number_density"=>Float64(density),
        "target_average_nup_per_cell"=>numerator(number),"individual_physical_sites_per_cell"=>n,
        "site_qn_up"=>up,"site_qn_down"=>down,"site_qn_name"=>"Sz",
        "physical_sz_background_per_cell"=>numerator(number)-n/2,
        "interpretation"=>"A rational mean filling of the infinite U1 ansatz, not separate conservation of each cell's physical number. Site Sz matrices retain eigenvalues +/-1/2.")
end
