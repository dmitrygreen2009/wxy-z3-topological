# Exact closed form of B† (sum_a W_ai μ_a^-) B on all eight matter states.
# This is an alternative representation of the microscopic Hamiltonian.
function real_charge_basis_opsum(lat)
    os=OpSum();ip=invperm(lat.order)
    a=(1-inv(sqrt(3)))/2
    b=2*(1+inv(sqrt(3)))
    d=2/sqrt(3)
    for v=1:lat.nv,i=1:3
        j,k=filter(!=(i),collect(1:3))
        si=ip[3(v-1)+i];sj=ip[3(v-1)+j];sk=ip[3(v-1)+k]
        sg=ip[3lat.nv+lat.legs[v][i]]
        add!(os,-a,"S-",si,"S+",sg)
        add!(os,-a,"S+",si,"S-",sg)
        add!(os,-b,"S-",si,"Sz",sj,"Sz",sk,"S+",sg)
        add!(os,-b,"S+",si,"Sz",sj,"Sz",sk,"S-",sg)
        add!(os,-d,"S+",si,"S-",sj,"S-",sk,"S+",sg)
        add!(os,-d,"S-",si,"S+",sj,"S+",sk,"S-",sg)
    end
    os
end
