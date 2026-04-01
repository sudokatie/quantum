@testset "Circuit Execution" begin
    @testset "empty circuit" begin
        c = Circuit(2)
        state = execute(c)
        @test state == zero_state(2)
    end
    
    @testset "single gate" begin
        # X on |0> gives |1>
        c = Circuit(1)
        add!(c, X, 1)
        state = execute(c)
        @test isapprox(state[1], 0, atol=1e-10)
        @test isapprox(state[2], 1, atol=1e-10)
        
        # H on |0> gives |+>
        c = Circuit(1)
        add!(c, H, 1)
        state = execute(c)
        @test isapprox(state[1], 1/sqrt(2), atol=1e-10)
        @test isapprox(state[2], 1/sqrt(2), atol=1e-10)
    end
    
    @testset "multi gate" begin
        # HH = I
        c = Circuit(1)
        add!(c, H, 1)
        add!(c, H, 1)
        state = execute(c)
        @test isapprox(state[1], 1, atol=1e-10)
        @test isapprox(state[2], 0, atol=1e-10)
        
        # XX = I
        c = Circuit(1)
        add!(c, X, 1)
        add!(c, X, 1)
        state = execute(c)
        @test isapprox(state[1], 1, atol=1e-10)
    end
    
    @testset "Bell state circuit" begin
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        
        state = execute(c)
        
        # Should be (|00> + |11>)/sqrt(2)
        @test isapprox(abs(state[1]), 1/sqrt(2), atol=1e-10)  # |00>
        @test isapprox(abs(state[2]), 0, atol=1e-10)          # |01>
        @test isapprox(abs(state[3]), 0, atol=1e-10)          # |10>
        @test isapprox(abs(state[4]), 1/sqrt(2), atol=1e-10)  # |11>
        @test is_normalized(state)
    end
    
    @testset "GHZ state" begin
        # |000> + |111> via H-CNOT-CNOT
        c = Circuit(3)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        add!(c, CNOT, 1, 3)
        
        state = execute(c)
        
        @test isapprox(abs(state[1]), 1/sqrt(2), atol=1e-10)  # |000>
        @test isapprox(abs(state[8]), 1/sqrt(2), atol=1e-10)  # |111>
        @test is_normalized(state)
    end
    
    @testset "custom initial state" begin
        c = Circuit(1)
        add!(c, X, 1)
        
        # Start from |1>
        state = execute(c, one_state(1))
        @test isapprox(state[1], 1, atol=1e-10)  # X|1> = |0>
    end
    
    @testset "controlled gate execution" begin
        c = Circuit(2)
        add_controlled!(c, X, [1], [2])
        
        # |00> unchanged
        state = execute(c, basis_state(2, 0))
        @test isapprox(state[1], 1, atol=1e-10)
        
        # |10> -> |11>
        state = execute(c, basis_state(2, 2))
        @test isapprox(state[4], 1, atol=1e-10)
    end
    
    @testset "execute_and_measure" begin
        # Deterministic: |0> circuit
        c = Circuit(1)
        counts = execute_and_measure(c; shots=100)
        @test counts[0] == 100
        
        # Probabilistic: H gate gives ~50/50
        c = Circuit(1)
        add!(c, H, 1)
        counts = execute_and_measure(c; shots=1000)
        @test haskey(counts, 0)
        @test haskey(counts, 1)
        @test 400 <= counts[0] <= 600
    end
    
    @testset "simulate" begin
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        
        counts = simulate(c, 1000)
        
        # Bell state: only |00> and |11>
        @test !haskey(counts, 1)  # no |01>
        @test !haskey(counts, 2)  # no |10>
        @test haskey(counts, 0)   # |00>
        @test haskey(counts, 3)   # |11>
    end
end
