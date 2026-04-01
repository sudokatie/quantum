@testset "Visualization" begin
    @testset "show_circuit" begin
        # Empty circuit
        c = Circuit(2)
        buf = IOBuffer()
        show_circuit(buf, c)
        output = String(take!(buf))
        @test contains(output, "Empty circuit")
        
        # Single gate circuit
        c = Circuit(2)
        add!(c, H, 1)
        buf = IOBuffer()
        show_circuit(buf, c)
        output = String(take!(buf))
        @test contains(output, "q1:")
        @test contains(output, "q2:")
        @test contains(output, "H")
        
        # Bell state circuit
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        buf = IOBuffer()
        show_circuit(buf, c)
        output = String(take!(buf))
        @test contains(output, "●")  # Control
        @test contains(output, "⊕")  # Target
    end
    
    @testset "show_state" begin
        # |0> state
        s = zero_state(1)
        buf = IOBuffer()
        show_state(buf, s)
        output = String(take!(buf))
        @test contains(output, "|0⟩")
        @test contains(output, "1.0")
        
        # |+> state
        s = apply(zero_state(1), H, 1)
        buf = IOBuffer()
        show_state(buf, s)
        output = String(take!(buf))
        @test contains(output, "|0⟩")
        @test contains(output, "|1⟩")
        
        # Bell state
        c = Circuit(2)
        add!(c, H, 1)
        add!(c, CNOT, 1, 2)
        s = execute(c)
        buf = IOBuffer()
        show_state(buf, s)
        output = String(take!(buf))
        @test contains(output, "|00⟩")
        @test contains(output, "|11⟩")
        @test !contains(output, "|01⟩")
        @test !contains(output, "|10⟩")
    end
    
    @testset "show_probabilities" begin
        # |0> has 100% probability of 0
        s = zero_state(1)
        buf = IOBuffer()
        show_probabilities(buf, s)
        output = String(take!(buf))
        @test contains(output, "|0⟩")
        @test contains(output, "100")
        
        # |+> has 50% each
        s = apply(zero_state(1), H, 1)
        buf = IOBuffer()
        show_probabilities(buf, s)
        output = String(take!(buf))
        @test contains(output, "50")
        @test contains(output, "█")
    end
    
    @testset "show_bloch" begin
        # |0> is at north pole (z=1)
        s = zero_state(1)
        buf = IOBuffer()
        show_bloch(buf, s)
        output = String(take!(buf))
        @test contains(output, "z = 1.0")
        
        # |1> is at south pole (z=-1)
        s = one_state(1)
        buf = IOBuffer()
        show_bloch(buf, s)
        output = String(take!(buf))
        @test contains(output, "z = -1.0")
        
        # |+> is on equator (z=0)
        s = apply(zero_state(1), H, 1)
        buf = IOBuffer()
        show_bloch(buf, s)
        output = String(take!(buf))
        @test contains(output, "z = 0.0")
        @test contains(output, "x = 1.0")
    end
end
