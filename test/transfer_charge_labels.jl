using Test,ITensors,KrylovKit
include("../src/transfer_charge_labels.jl")
@testset "Transfer vector charge labels in the original QN basis" begin
    i=Index([QN("Sz",-1)=>1,QN("Sz",1)=>1],"virtual";dir=ITensors.Out)
    original=[i,dag(prime(i))];plain=ITensors.removeqns.(original)
    for (a,b,charge) in [(1,1,0),(2,2,0),(1,2,-2),(2,1,2)]
        matrix=zeros(ComplexF64,2,2);matrix[a,b]=1
        weights=transfer_charge_weights(ITensor(matrix,plain...),original)
        @test weights[charge]≈1 atol=1e-15
        @test sum(values(weights))≈1 atol=1e-15
    end
    matrix=ComplexF64[0 1;im 0]/sqrt(2)
    weights=transfer_charge_weights(ITensor(matrix,plain...),original)
    @test weights[-2]≈.5 atol=1e-15
    @test weights[2]≈.5 atol=1e-15
end
