function local_eigenpair_record(M,lambda,v,info,tol)
    residual=norm(M(v)-lambda*v)/norm(v)
    reported=info.normres isa Number ? abs(info.normres) : maximum(abs.(info.normres);init=0.)
    Dict("operator_type"=>string(typeof(M)),"operator_site"=>hasproperty(M,:n) ? getproperty(M,:n) : nothing,
        "eigenvalue_real"=>real(lambda),"eigenvalue_imag"=>imag(lambda),
        "requested_tolerance"=>tol,"true_eigenpair_residual"=>residual,
        "scaled_eigenpair_residual"=>residual/max(1.,abs(lambda)),
        "reported_normres_max"=>reported,"converged_eigenpairs"=>info.converged,
        "krylov_iterations"=>info.numiter,"operator_applications"=>info.numops,
        "local_criterion_passed"=>info.converged>=1 && residual<=tol)
end
