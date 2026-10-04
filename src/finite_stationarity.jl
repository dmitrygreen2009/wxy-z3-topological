"""Post-update/saved-state two-site projected residual using library ProjMPO.
This includes truncation and variational-basis limitations; it is not the
pre-truncation Krylov residual or a certificate of the global ground state.
"""
function finite_bond_stationarity(H,psi,b)
    ITensorMPS.orthogonalize!(psi,b)
    P=ProjMPO(H);ITensorMPS.position!(P,psi,b)
    v=psi[b]*psi[b+1];Mv=P(v);q=inner(v,Mv)/inner(v,v)
    Dict("bond"=>b,"projected_rayleigh_energy_real"=>real(q),
        "projected_rayleigh_energy_imag"=>imag(q),
        "projected_residual"=>norm(Mv-q*v)/norm(v),
        "interpretation"=>"Residual in the retained two-site variational basis; vanishing residual does not certify the full Hilbert-space eigenstate or ground state.")
end
