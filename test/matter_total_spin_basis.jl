using Test
isdefined(Main,:cylinder) || include("../src/model.jl")
include("../src/matter_total_spin_basis.jl")
@testset "Exact quartet and cycle-multiplicity doublets use the validated MPS encoder" begin
    data=matter_total_spin_basis();U=data.U
    independent=JSON3.read(read("geometry/matter_total_spin_basis.json",String),Dict{String,Any})
    reference=[complex(independent["unitary_real"][i][j],independent["unitary_imag"][i][j]) for i=1:8,j=1:8]
    @test U≈reference atol=1e-13
    @test U'U≈Matrix{ComplexF64}(I,8,8) atol=1e-13
    C=zeros(ComplexF64,8,8);Sp=zeros(ComplexF64,8,8)
    for old=0:7
        new=sum(((old>>a)&1)<<mod(a+1,3) for a=0:2);C[new+1,old+1]=1
        for a=0:2
            iszero((old>>a)&1) && (Sp[(old|(1<<a))+1,old+1]=1)
        end
    end
    Sz=Diagonal([count_ones(b)-1.5 for b=0:7]);S2=Sz^2+(Sp*Sp'+Sp'*Sp)/2
    @test U'S2*U≈Diagonal([label.S*(label.S+1) for label in data.labels]) atol=1e-13
    @test U'Sz*U≈Diagonal([label.m for label in data.labels]) atol=1e-13
    @test U'C*U≈Diagonal([cis(2pi*label.k/3) for label in data.labels]) atol=1e-13
    for (b,label) in zip(data.encoded_bits,data.labels)
        @test count_ones(b)==label.m+1.5
        @test mod(sum(a*((b>>a)&1) for a=0:2),3)==label.k
    end
end
