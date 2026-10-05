# Exact 4⊕2⊕2 decomposition, expressed through the already validated encoder.
isdefined(@__MODULE__,:matter_charge_basis) || include("matter_charge_basis.jl")
function matter_total_spin_basis()
    encoded_bits=[0,1,6,7,2,3,4,5]
    phases=[1,1,1,1,1,-1,1,-1]
    U=matter_charge_basis()[:,encoded_bits.+1]*Diagonal(phases)
    labels=[(;S=3/2,m=-3/2,k=0),(;S=3/2,m=-1/2,k=0),
        (;S=3/2,m=1/2,k=0),(;S=3/2,m=3/2,k=0),
        (;S=1/2,m=-1/2,k=1),(;S=1/2,m=1/2,k=1),
        (;S=1/2,m=-1/2,k=2),(;S=1/2,m=1/2,k=2)]
    (;U,labels,encoded_bits,phases)
end
