@testset "Circuit" begin
    @testset "construction" begin
        c = Circuit(3)
        @test qubits(c) == 3
        @test gate_count(c) == 0
        @test depth(c) == 0
    end
    
    @testset "add single qubit gates" begin
        c = Circuit(2)
        add!(c, H, 1)
        @test gate_count(c) == 1
        @test c.operations[1].name == "H"
        @test c.operations[1].targets == [1]
        
        add!(c, X, 2)
        @test gate_count(c) == 2
        @test c.operations[2].name == "X"
    end
    
    @testset "add two qubit gates" begin
        c = Circuit(2)
        add!(c, CNOT, 1, 2)
        @test gate_count(c) == 1
        @test c.operations[1].name == "CNOT"
        @test c.operations[1].targets == [1, 2]
        
        add!(c, SWAP, 1, 2)
        @test gate_count(c) == 2
        @test c.operations[2].name == "SWAP"
    end
    
    @testset "add_controlled!" begin
        c = Circuit(3)
        add_controlled!(c, X, [1], [2])
        @test gate_count(c) == 1
        @test c.operations[1].controls == [1]
        @test c.operations[1].targets == [2]
        @test c.operations[1].name == "CX"
        
        # Multiple controls
        add_controlled!(c, X, [1, 2], [3])
        @test c.operations[2].controls == [1, 2]
        @test c.operations[2].targets == [3]
    end
    
    @testset "depth calculation" begin
        # Single gate: depth 1
        c = Circuit(2)
        add!(c, H, 1)
        @test depth(c) == 1
        
        # Parallel gates: still depth 1
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, H, 2)
        @test depth(c) == 1
        
        # Sequential gates on same qubit: depth 2
        c = Circuit(1)
        add!(c, H, 1)
        add!(c, X, 1)
        @test depth(c) == 2
        
        # CNOT adds depth
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        @test depth(c) == 2
    end
    
    @testset "gate_count" begin
        c = Circuit(3)
        @test gate_count(c) == 0
        
        add!(c, H, 1)
        add!(c, H, 2)
        add!(c, H, 3)
        @test gate_count(c) == 3
        
        add!(c, CNOT, 1, 2)
        @test gate_count(c) == 4
    end
    
    @testset "convenience methods" begin
        c = Circuit(2)
        add_h!(c, 1)
        add_x!(c, 2)
        add_cnot!(c, 1, 2)
        
        @test gate_count(c) == 3
        @test c.operations[1].name == "H"
        @test c.operations[2].name == "X"
        @test c.operations[3].name == "CNOT"
    end
    
    @testset "bell state circuit" begin
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        
        @test gate_count(c) == 2
        @test depth(c) == 2
    end
    
    @testset "out of range errors" begin
        c = Circuit(2)
        @test_throws AssertionError add!(c, H, 3)
        @test_throws AssertionError add!(c, H, 0)
        @test_throws AssertionError add_controlled!(c, X, [3], [1])
    end
end
