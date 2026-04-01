using Test
using Quantum
using LinearAlgebra: I

@testset "Quantum.jl" begin
    include("test_state.jl")
    include("test_gates.jl")
    include("test_apply.jl")
    include("test_measure.jl")
end
