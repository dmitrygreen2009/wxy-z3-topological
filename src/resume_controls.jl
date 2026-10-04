# Preserve explicitly recorded local DMRG controls on checkpoint continuation.
function finite_resume_controls(settings,phase)
    dimension=occursin("refinement",phase) ? get(settings,"refinement_krylovdim",get(settings,"krylovdim",12)) : get(settings,"krylovdim",12)
    (;eigsolve_krylovdim=dimension,eigsolve_maxiter=get(settings,"eigsolve_maxiter",1),
      eigsolve_tol=get(settings,"eigsolve_tol",1e-14),eigsolve_verbosity=get(settings,"eigsolve_verbosity",0))
end
