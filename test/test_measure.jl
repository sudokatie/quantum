@testset "Measurement" begin
    @testset "probabilities" begin
        # |0> has P(0)=1, P(1)=0
        s = zero_state(1)
        p = probabilities(s)
        @test isapprox(p[1], 1.0, atol=1e-10)
        @test isapprox(p[2], 0.0, atol=1e-10)
        
        # |1> has P(0)=0, P(1)=1
        s = one_state(1)
        p = probabilities(s)
        @test isapprox(p[1], 0.0, atol=1e-10)
        @test isapprox(p[2], 1.0, atol=1e-10)
        
        # |+> has P(0)=0.5, P(1)=0.5
        s = apply(zero_state(1), H, 1)
        p = probabilities(s)
        @test isapprox(p[1], 0.5, atol=1e-10)
        @test isapprox(p[2], 0.5, atol=1e-10)
        
        # Probabilities sum to 1
        s = zero_state(3)
        s = apply(s, H, 1)
        s = apply(s, H, 2)
        p = probabilities(s)
        @test isapprox(sum(p), 1.0, atol=1e-10)
    end
    
    @testset "measure deterministic" begin
        # |0> always measures 0
        s = zero_state(1)
        for _ in 1:10
            result, collapsed = measure(s)
            @test result == 0
            @test collapsed == zero_state(1)
        end
        
        # |1> always measures 1
        s = one_state(1)
        for _ in 1:10
            result, collapsed = measure(s)
            @test result == 1
            @test collapsed == one_state(1)
        end
        
        # |11> always measures 3
        s = basis_state(2, 3)
        for _ in 1:10
            result, _ = measure(s)
            @test result == 3
        end
    end
    
    @testset "measure_qubit" begin
        # Measure qubit 1 of |0> gives 0
        s = zero_state(1)
        bit, collapsed = measure_qubit(s, 1)
        @test bit == 0
        @test is_normalized(collapsed)
        
        # Measure qubit 1 of |1> gives 1
        s = one_state(1)
        bit, collapsed = measure_qubit(s, 1)
        @test bit == 1
        @test is_normalized(collapsed)
        
        # Measure qubit 1 of |10> gives 1
        s = basis_state(2, 2)  # |10>
        bit, collapsed = measure_qubit(s, 1)
        @test bit == 1
        @test is_normalized(collapsed)
        
        # Measure qubit 2 of |10> gives 0
        s = basis_state(2, 2)  # |10>
        bit, collapsed = measure_qubit(s, 2)
        @test bit == 0
        @test is_normalized(collapsed)
    end
    
    @testset "measure_qubit superposition" begin
        # Bell state: measuring qubit 1 collapses both
        # |Phi+> = (|00> + |11>)/sqrt(2)
        s = zero_state(2)
        s = apply(s, H, 1)
        s = apply(s, CNOT, [1, 2])
        
        bit, collapsed = measure_qubit(s, 1)
        @test bit in [0, 1]
        @test is_normalized(collapsed)
        
        # After measuring qubit 1, the state is collapsed
        if bit == 0
            @test isapprox(abs(collapsed[1]), 1, atol=1e-10)  # |00>
        else
            @test isapprox(abs(collapsed[4]), 1, atol=1e-10)  # |11>
        end
    end
    
    @testset "sample" begin
        # Sample |0> should always give 0
        s = zero_state(1)
        counts = sample(s, 100)
        @test length(counts) == 1
        @test counts[0] == 100
        
        # Sample |+> should give roughly 50/50
        s = apply(zero_state(1), H, 1)
        counts = sample(s, 1000)
        @test haskey(counts, 0)
        @test haskey(counts, 1)
        # Allow for statistical fluctuation
        @test 400 <= get(counts, 0, 0) <= 600
        @test 400 <= get(counts, 1, 0) <= 600
    end
    
    @testset "expectation" begin
        # <0|Z|0> = 1
        s = zero_state(1)
        @test isapprox(expectation(s, Z), 1.0, atol=1e-10)
        
        # <1|Z|1> = -1
        s = one_state(1)
        @test isapprox(expectation(s, Z), -1.0, atol=1e-10)
        
        # <+|Z|+> = 0
        s = apply(zero_state(1), H, 1)
        @test isapprox(expectation(s, Z), 0.0, atol=1e-10)
        
        # <+|X|+> = 1
        s = apply(zero_state(1), H, 1)
        @test isapprox(expectation(s, X), 1.0, atol=1e-10)
        
        # <0|X|0> = 0
        s = zero_state(1)
        @test isapprox(expectation(s, X), 0.0, atol=1e-10)
    end
end
