@testset "Examples" begin
    @testset "bell_state example runs" begin
        # Verify Bell state creation works
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        state = execute(c)
        
        # Bell state should have concurrence 1
        @test isapprox(concurrence(state), 1.0, atol=1e-6)
        
        # Only |00⟩ and |11⟩ have probability
        probs = probabilities(state)
        @test isapprox(probs[1], 0.5, atol=1e-10)  # |00⟩
        @test isapprox(probs[2], 0.0, atol=1e-10)  # |01⟩
        @test isapprox(probs[3], 0.0, atol=1e-10)  # |10⟩
        @test isapprox(probs[4], 0.5, atol=1e-10)  # |11⟩
    end
    
    @testset "grover example runs" begin
        # Simplified 2-qubit Grover for target |11⟩
        c = Circuit(2)
        
        # Superposition
        add!(c, H, 1)
        add!(c, H, 2)
        
        # Oracle (CZ marks |11⟩)
        add!(c, CZ, 1, 2)
        
        # Diffusion
        add!(c, H, 1)
        add!(c, H, 2)
        add!(c, X, 1)
        add!(c, X, 2)
        add!(c, CZ, 1, 2)
        add!(c, X, 1)
        add!(c, X, 2)
        add!(c, H, 1)
        add!(c, H, 2)
        
        state = execute(c)
        probs = probabilities(state)
        
        # Target |11⟩ should have highest probability
        @test probs[4] > probs[1]  # |11⟩ > |00⟩
        @test probs[4] > probs[2]  # |11⟩ > |01⟩
        @test probs[4] > probs[3]  # |11⟩ > |10⟩
    end
    
    @testset "qft example runs" begin
        # QFT on |000⟩ should give uniform superposition
        n = 3
        c = Circuit(n)
        
        # Simple QFT without swaps for testing
        for j in 1:n
            add!(c, H, j)
        end
        
        state = execute(c)
        probs = probabilities(state)
        
        # All 8 states should have equal probability
        expected = 1.0 / 2^n
        for p in probs
            @test isapprox(p, expected, atol=1e-10)
        end
    end
    
    @testset "teleportation setup" begin
        # Test that Bell pair is created correctly
        c = Circuit(3)
        add!(c, H, 2)
        add!(c, CNOT, 2, 3)
        
        initial = tensor(zero_state(1), zero_state(2))
        state = execute(c, initial)
        
        # Should be |0⟩|Φ+⟩
        @test is_normalized(state)
        @test length(state) == 8
    end
end
