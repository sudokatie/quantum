using Test
using Quantum
using LinearAlgebra: I

@testset "Quantum.jl" begin
    include("test_state.jl")
    include("test_gates.jl")
    include("test_apply.jl")
    include("test_measure.jl")
    include("test_entangle.jl")
    include("test_circuit.jl")
    include("test_execute.jl")
    include("test_visualize.jl")
    include("test_examples.jl")
    include("test_optimizer.jl")
    include("test_qpe.jl")
end
