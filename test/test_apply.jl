@testset "Gate Application" begin
    @testset "Single qubit gates" begin
        # X|0> = |1>
        s = zero_state(1)
        s = apply(s, X, 1)
        @test isapprox(s[1], 0, atol=1e-10)
        @test isapprox(s[2], 1, atol=1e-10)
        
        # H|0> = |+>
        s = zero_state(1)
        s = apply(s, H, 1)
        @test isapprox(s[1], 1/sqrt(2), atol=1e-10)
        @test isapprox(s[2], 1/sqrt(2), atol=1e-10)
        
        # Z|1> = -|1>
        s = one_state(1)
        s = apply(s, Z, 1)
        @test isapprox(s[1], 0, atol=1e-10)
        @test isapprox(s[2], -1, atol=1e-10)
    end
    
    @testset "Multi-qubit state single gate" begin
        # Apply X to qubit 1 of |00>
        s = zero_state(2)
        s = apply(s, X, 1)
        # Result should be |10>
        @test isapprox(s[1], 0, atol=1e-10)  # |00>
        @test isapprox(s[2], 0, atol=1e-10)  # |01>
        @test isapprox(s[3], 1, atol=1e-10)  # |10>
        @test isapprox(s[4], 0, atol=1e-10)  # |11>
        
        # Apply X to qubit 2 of |00>
        s = zero_state(2)
        s = apply(s, X, 2)
        # Result should be |01>
        @test isapprox(s[1], 0, atol=1e-10)  # |00>
        @test isapprox(s[2], 1, atol=1e-10)  # |01>
        @test isapprox(s[3], 0, atol=1e-10)  # |10>
        @test isapprox(s[4], 0, atol=1e-10)  # |11>
    end
    
    @testset "Normalization preserved" begin
        s = zero_state(3)
        s = apply(s, H, 1)
        s = apply(s, H, 2)
        s = apply(s, H, 3)
        @test is_normalized(s)
    end
    
    @testset "Sequential gates" begin
        # HH = I
        s = zero_state(1)
        s = apply(s, H, 1)
        s = apply(s, H, 1)
        @test isapprox(s[1], 1, atol=1e-10)
        @test isapprox(s[2], 0, atol=1e-10)
        
        # XX = I
        s = zero_state(1)
        s = apply(s, X, 1)
        s = apply(s, X, 1)
        @test isapprox(s[1], 1, atol=1e-10)
        @test isapprox(s[2], 0, atol=1e-10)
    end
    
    @testset "CNOT application" begin
        # CNOT|00> = |00>
        s = zero_state(2)
        s = apply(s, CNOT, [1, 2])
        @test isapprox(s[1], 1, atol=1e-10)
        
        # CNOT|10> = |11>
        s = basis_state(2, 2)  # |10>
        s = apply(s, CNOT, [1, 2])
        @test isapprox(s[4], 1, atol=1e-10)  # |11>
        
        # CNOT|11> = |10>
        s = basis_state(2, 3)  # |11>
        s = apply(s, CNOT, [1, 2])
        @test isapprox(s[3], 1, atol=1e-10)  # |10>
    end
    
    @testset "Bell state creation" begin
        # |00> -> H on qubit 1 -> CNOT -> Bell state
        s = zero_state(2)
        s = apply(s, H, 1)
        s = apply(s, CNOT, [1, 2])
        
        # Should be (|00> + |11>)/sqrt(2)
        @test isapprox(abs(s[1]), 1/sqrt(2), atol=1e-10)  # |00>
        @test isapprox(abs(s[2]), 0, atol=1e-10)          # |01>
        @test isapprox(abs(s[3]), 0, atol=1e-10)          # |10>
        @test isapprox(abs(s[4]), 1/sqrt(2), atol=1e-10)  # |11>
        @test is_normalized(s)
    end
    
    @testset "In-place application" begin
        s = zero_state(2)
        apply!(s, X, 1)
        @test isapprox(s[3], 1, atol=1e-10)  # |10>
        
        apply!(s, H, 1)
        @test is_normalized(s)
    end
    
    @testset "Three qubit state" begin
        s = zero_state(3)
        s = apply(s, X, 2)  # |010>
        @test isapprox(s[3], 1, atol=1e-10)  # Index 3 = binary 010
        
        s = apply(s, X, 3)  # |011>
        @test isapprox(s[4], 1, atol=1e-10)  # Index 4 = binary 011
    end
    
    @testset "apply_controlled" begin
        # apply_controlled with single control = CNOT
        s = basis_state(2, 2)  # |10>
        s = apply_controlled(s, X, [1], [2])
        @test isapprox(s[4], 1, atol=1e-10)  # |11>
        
        # Control not set, no change
        s = basis_state(2, 0)  # |00>
        s = apply_controlled(s, X, [1], [2])
        @test isapprox(s[1], 1, atol=1e-10)  # Still |00>
        
        # Toffoli-like: two controls
        s = basis_state(3, 6)  # |110>
        s = apply_controlled(s, X, [1, 2], [3])
        @test isapprox(s[8], 1, atol=1e-10)  # |111>
        
        # Only one control set, no flip
        s = basis_state(3, 4)  # |100>
        s = apply_controlled(s, X, [1, 2], [3])
        @test isapprox(s[5], 1, atol=1e-10)  # Still |100>
    end
    
    @testset "qubit limit" begin
        # 20 qubits should work
        @test_nowarn zero_state(20)
        
        # 21 qubits should fail
        @test_throws ErrorException zero_state(21)
        @test_throws ErrorException one_state(21)
        @test_throws ErrorException basis_state(21, 0)
    end
end
