using LinearAlgebra: tr, eigvals

@testset "Entanglement" begin
    @testset "tensor product" begin
        # |0> ⊗ |0> = |00>
        a = zero_state(1)
        b = zero_state(1)
        ab = tensor(a, b)
        @test ab.n_qubits == 2
        @test isapprox(ab[1], 1, atol=1e-10)  # |00>
        
        # |0> ⊗ |1> = |01>
        a = zero_state(1)
        b = one_state(1)
        ab = tensor(a, b)
        @test isapprox(ab[2], 1, atol=1e-10)  # |01>
        
        # |1> ⊗ |0> = |10>
        a = one_state(1)
        b = zero_state(1)
        ab = tensor(a, b)
        @test isapprox(ab[3], 1, atol=1e-10)  # |10>
        
        # |+> ⊗ |0> = (|00> + |10>)/sqrt(2)
        plus = apply(zero_state(1), H, 1)
        ab = tensor(plus, zero_state(1))
        @test ab.n_qubits == 2
        @test isapprox(abs(ab[1]), 1/sqrt(2), atol=1e-10)
        @test isapprox(abs(ab[3]), 1/sqrt(2), atol=1e-10)
        @test is_normalized(ab)
    end
    
    @testset "partial_trace" begin
        # Product state |00>: trace out qubit 2
        s = zero_state(2)
        rho = partial_trace(s, [1])
        @test size(rho) == (2, 2)
        @test isapprox(rho[1,1], 1, atol=1e-10)  # |0><0|
        @test isapprox(rho[2,2], 0, atol=1e-10)
        
        # Product state |00>: trace out qubit 1
        rho = partial_trace(s, [2])
        @test isapprox(rho[1,1], 1, atol=1e-10)  # |0><0|
        
        # Bell state: trace out qubit 2 gives maximally mixed
        s = zero_state(2)
        s = apply(s, H, 1)
        s = apply(s, CNOT, [1, 2])
        
        rho = partial_trace(s, [1])
        @test isapprox(rho[1,1], 0.5, atol=1e-10)  # Maximally mixed
        @test isapprox(rho[2,2], 0.5, atol=1e-10)
        @test isapprox(abs(rho[1,2]), 0, atol=1e-10)
    end
    
    @testset "purity" begin
        # Pure state has purity 1
        s = zero_state(1)
        rho = s.amplitudes * s.amplitudes'
        @test isapprox(purity(rho), 1.0, atol=1e-10)
        
        # Maximally mixed 1-qubit state has purity 0.5
        mixed = ComplexF64[0.5 0; 0 0.5]
        @test isapprox(purity(mixed), 0.5, atol=1e-10)
    end
    
    @testset "concurrence" begin
        # Product state |00> has concurrence 0
        s = zero_state(2)
        @test isapprox(concurrence(s), 0.0, atol=1e-6)
        
        # Product state |+0> has concurrence 0
        plus = apply(zero_state(1), H, 1)
        s = tensor(plus, zero_state(1))
        @test isapprox(concurrence(s), 0.0, atol=1e-6)
        
        # Bell state has concurrence 1
        s = zero_state(2)
        s = apply(s, H, 1)
        s = apply(s, CNOT, [1, 2])
        @test isapprox(concurrence(s), 1.0, atol=1e-6)
    end
    
    @testset "is_separable" begin
        # Product states are separable
        s = zero_state(2)
        @test is_separable(s)
        
        s = basis_state(2, 3)  # |11>
        @test is_separable(s)
        
        # |+-> is separable (product of |+> and |->)
        plus = apply(zero_state(1), H, 1)
        minus = apply(one_state(1), H, 1)
        s = tensor(plus, minus)
        @test is_separable(s)
        
        # Bell state is not separable
        s = zero_state(2)
        s = apply(s, H, 1)
        s = apply(s, CNOT, [1, 2])
        @test !is_separable(s)
    end
    
    @testset "is_entangled" begin
        # Product states are not entangled
        s = zero_state(2)
        @test !is_entangled(s)
        
        s = basis_state(2, 3)  # |11>
        @test !is_entangled(s)
        
        # Bell state is entangled
        s = zero_state(2)
        s = apply(s, H, 1)
        s = apply(s, CNOT, [1, 2])
        @test is_entangled(s)
        
        # is_entangled should be opposite of is_separable
        plus = apply(zero_state(1), H, 1)
        s = tensor(plus, zero_state(1))
        @test is_separable(s) == !is_entangled(s)
    end
    
    @testset "von_neumann_entropy" begin
        # Pure state has entropy 0
        rho = ComplexF64[1 0; 0 0]
        @test isapprox(von_neumann_entropy(rho), 0.0, atol=1e-10)
        
        # Maximally mixed 1-qubit state has entropy 1
        mixed = ComplexF64[0.5 0; 0 0.5]
        @test isapprox(von_neumann_entropy(mixed), 1.0, atol=1e-10)
        
        # Bell state reduced to 1 qubit is maximally mixed (entropy 1)
        s = zero_state(2)
        s = apply(s, H, 1)
        s = apply(s, CNOT, [1, 2])
        rho = partial_trace(s, [1])
        @test isapprox(von_neumann_entropy(rho), 1.0, atol=1e-6)
    end
end
