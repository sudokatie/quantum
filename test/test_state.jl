@testset "State Vector" begin
    @testset "zero_state" begin
        s = zero_state(1)
        @test s.n_qubits == 1
        @test s[1] == 1.0
        @test s[2] == 0.0
        
        s2 = zero_state(2)
        @test s2.n_qubits == 2
        @test length(s2) == 4
        @test s2[1] == 1.0
    end
    
    @testset "one_state" begin
        s = one_state(1)
        @test s[1] == 0.0
        @test s[2] == 1.0
        
        s2 = one_state(2)
        @test s2[4] == 1.0
    end
    
    @testset "basis_state" begin
        s = basis_state(2, 0)  # |00>
        @test s[1] == 1.0
        
        s = basis_state(2, 3)  # |11>
        @test s[4] == 1.0
        
        @test_throws ErrorException basis_state(2, 4)
        @test_throws ErrorException basis_state(2, -1)
    end
    
    @testset "normalization" begin
        s = zero_state(1)
        @test is_normalized(s)
        
        unnorm = StateVector(ComplexF64[2.0, 0.0])
        @test !is_normalized(unnorm)
        normalize!(unnorm)
        @test is_normalized(unnorm)
        @test abs(unnorm[1] - 1.0) < 1e-10
    end
    
    @testset "copy and equality" begin
        s = zero_state(2)
        s2 = copy(s)
        @test s == s2
        @test s !== s2  # Different objects
    end
    
    @testset "num_qubits" begin
        @test num_qubits(zero_state(1)) == 1
        @test num_qubits(zero_state(3)) == 3
        @test num_qubits(zero_state(5)) == 5
    end
    
    @testset "invalid construction" begin
        @test_throws ErrorException StateVector(ComplexF64[1.0, 0.0, 0.0])
        @test_throws ErrorException StateVector(ComplexF64[1.0, 0.0, 0.0, 0.0, 0.0])
    end
end
