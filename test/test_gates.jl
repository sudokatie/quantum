@testset "Quantum Gates" begin
    @testset "Gate dimensions" begin
        @test size(I_GATE) == (2, 2)
        @test size(X) == (2, 2)
        @test size(Y) == (2, 2)
        @test size(Z) == (2, 2)
        @test size(H) == (2, 2)
        @test size(S) == (2, 2)
        @test size(T) == (2, 2)
    end
    
    @testset "Unitarity" begin
        @test is_unitary(I_GATE)
        @test is_unitary(X)
        @test is_unitary(Y)
        @test is_unitary(Z)
        @test is_unitary(H)
        @test is_unitary(S)
        @test is_unitary(T)
    end
    
    @testset "Pauli gates" begin
        # X|0> = |1>
        zero = ComplexF64[1, 0]
        one = ComplexF64[0, 1]
        @test X * zero == one
        @test X * one == zero
        
        # Z|0> = |0>, Z|1> = -|1>
        @test Z * zero == zero
        @test Z * one == -one
        
        # X^2 = Y^2 = Z^2 = I
        @test isapprox(X * X, I_GATE, atol=1e-10)
        @test isapprox(Y * Y, I_GATE, atol=1e-10)
        @test isapprox(Z * Z, I_GATE, atol=1e-10)
    end
    
    @testset "Hadamard" begin
        zero = ComplexF64[1, 0]
        one = ComplexF64[0, 1]
        
        # H|0> = |+> = (|0> + |1>)/sqrt(2)
        plus = H * zero
        @test isapprox(plus[1], 1/sqrt(2), atol=1e-10)
        @test isapprox(plus[2], 1/sqrt(2), atol=1e-10)
        
        # H|1> = |-> = (|0> - |1>)/sqrt(2)
        minus = H * one
        @test isapprox(minus[1], 1/sqrt(2), atol=1e-10)
        @test isapprox(minus[2], -1/sqrt(2), atol=1e-10)
        
        # H^2 = I
        @test isapprox(H * H, I_GATE, atol=1e-10)
    end
    
    @testset "Rotation gates" begin
        # Rx(0) = Ry(0) = Rz(0) = I (up to global phase)
        @test isapprox(Rx(0), I_GATE, atol=1e-10)
        @test isapprox(Ry(0), I_GATE, atol=1e-10)
        # Rz(0) gives [exp(0), exp(0)] = I
        @test isapprox(Rz(0), I_GATE, atol=1e-10)
        
        # Unitarity of rotation gates
        @test is_unitary(Rx(pi/4))
        @test is_unitary(Ry(pi/3))
        @test is_unitary(Rz(pi/6))
        
        # Rx(pi) = -i*X (up to global phase, check matrix structure)
        rx_pi = Rx(pi)
        @test isapprox(abs(rx_pi[1,1]), 0, atol=1e-10)
        @test isapprox(abs(rx_pi[1,2]), 1, atol=1e-10)
    end
    
    @testset "Two-qubit gates" begin
        @test size(CNOT) == (4, 4)
        @test size(CZ) == (4, 4)
        @test size(SWAP) == (4, 4)
        @test size(iSWAP) == (4, 4)
        
        @test is_unitary(CNOT)
        @test is_unitary(CZ)
        @test is_unitary(SWAP)
        @test is_unitary(iSWAP)
    end
    
    @testset "CNOT behavior" begin
        # |00> -> |00>
        ket00 = ComplexF64[1, 0, 0, 0]
        @test CNOT * ket00 == ket00
        
        # |01> -> |01>
        ket01 = ComplexF64[0, 1, 0, 0]
        @test CNOT * ket01 == ket01
        
        # |10> -> |11>
        ket10 = ComplexF64[0, 0, 1, 0]
        ket11 = ComplexF64[0, 0, 0, 1]
        @test CNOT * ket10 == ket11
        
        # |11> -> |10>
        @test CNOT * ket11 == ket10
    end
    
    @testset "SWAP behavior" begin
        ket01 = ComplexF64[0, 1, 0, 0]
        ket10 = ComplexF64[0, 0, 1, 0]
        
        @test SWAP * ket01 == ket10
        @test SWAP * ket10 == ket01
    end
    
    @testset "CZ behavior" begin
        ket11 = ComplexF64[0, 0, 0, 1]
        @test CZ * ket11 == -ket11
        
        # Other states unchanged
        ket00 = ComplexF64[1, 0, 0, 0]
        @test CZ * ket00 == ket00
    end
    
    @testset "controlled gate construction" begin
        # controlled(X) should equal CNOT
        cx = controlled(X)
        @test isapprox(cx, CNOT, atol=1e-10)
        
        # controlled(Z) should equal CZ
        cz = controlled(Z)
        @test isapprox(cz, CZ, atol=1e-10)
    end
    
    @testset "tensor_product" begin
        # I (x) I should be 4x4 identity
        ii = tensor_product(I_GATE, I_GATE)
        @test size(ii) == (4, 4)
        @test isapprox(ii, Matrix{ComplexF64}(I, 4, 4), atol=1e-10)
        
        # Dimension check
        hh = tensor_product(H, H)
        @test size(hh) == (4, 4)
        @test is_unitary(hh)
    end
end
