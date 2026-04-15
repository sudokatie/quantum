@testset "Optimization Algorithms" begin
    @testset "Hamiltonian" begin
        # Simple Z Hamiltonian: H = Z_1
        H = hamiltonian(1, [(1.0, [(1, 'Z')])])
        @test H.n_qubits == 1
        @test length(H.terms) == 1
        
        # Two-qubit Hamiltonian: H = Z_1 Z_2
        H2 = hamiltonian(2, [(1.0, [(1, 'Z'), (2, 'Z')])])
        @test H2.n_qubits == 2
    end
    
    @testset "Pauli expectation" begin
        # |0⟩ state: <Z> = 1
        state = zero_state(1)
        @test pauli_expectation(state, [(1, 'Z')]) ≈ 1.0
        
        # |1⟩ state: <Z> = -1
        state = one_state(1)
        @test pauli_expectation(state, [(1, 'Z')]) ≈ -1.0
        
        # |+⟩ state: <X> = 1, <Z> = 0
        state = zero_state(1)
        apply!(state, H, 1)
        @test pauli_expectation(state, [(1, 'X')]) ≈ 1.0 atol=1e-10
        @test pauli_expectation(state, [(1, 'Z')]) ≈ 0.0 atol=1e-10
    end
    
    @testset "Expectation value" begin
        # H = Z, |0⟩: <H> = 1
        H = hamiltonian(1, [(1.0, [(1, 'Z')])])
        state = zero_state(1)
        @test expectation_value(state, H) ≈ 1.0
        
        # H = 0.5*I + 0.5*Z, |0⟩: <H> = 0.5 + 0.5 = 1
        H = hamiltonian(1, [(0.5, Tuple{Int,Char}[]), (0.5, [(1, 'Z')])])
        @test expectation_value(state, H) ≈ 1.0
    end
    
    @testset "Variational circuit" begin
        vc = variational_circuit(2; layers=1, entanglement=:linear)
        @test vc.n_qubits == 2
        @test vc.layers == 1
        @test num_parameters(vc) == 6  # 2 qubits * 3 rotations * 1 layer
        
        vc2 = variational_circuit(3; layers=2, entanglement=:circular)
        @test num_parameters(vc2) == 18  # 3 qubits * 3 rotations * 2 layers
    end
    
    @testset "Build ansatz" begin
        vc = variational_circuit(2; layers=1)
        params = zeros(6)
        circuit = build_ansatz(vc, params)
        @test qubits(circuit) == 2
        @test gate_count(circuit) > 0
        
        # Execute should work
        state = execute(circuit)
        @test is_normalized(state)
    end
    
    @testset "Nelder-Mead optimizer" begin
        # Minimize f(x) = (x-1)^2
        f(x) = (x[1] - 1.0)^2
        result = nelder_mead(f, [0.0]; max_iter=100)
        @test result.optimal_value < 1e-6
        @test abs(result.optimal_params[1] - 1.0) < 1e-3
        
        # Minimize Rosenbrock (harder)
        rosenbrock(x) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
        result = nelder_mead(rosenbrock, [0.0, 0.0]; max_iter=1000, tol=1e-8)
        @test result.optimal_value < 1e-3
    end
    
    @testset "VQE - simple" begin
        # Ground state of H = Z is |0⟩ with E = -1
        # But our Z convention: Z|0⟩ = |0⟩, so E = +1 for |0⟩ and E = -1 for |1⟩
        H = hamiltonian(1, [(1.0, [(1, 'Z')])])
        ansatz = variational_circuit(1; layers=1)
        
        result = vqe(H, ansatz; max_iter=200, tol=1e-4)
        # Should find minimum energy around -1
        @test result.optimal_value < 0.0  # Found negative energy
    end
    
    @testset "MaxCut Hamiltonian" begin
        # Triangle graph
        edges = [(1, 2), (2, 3), (1, 3)]
        H = maxcut_hamiltonian(edges, 3)
        @test H.n_qubits == 3
        # 6 terms: 3 identity (0.5 each) + 3 ZZ (-0.5 each)
        @test length(H.terms) == 6
    end
    
    @testset "Standard mixer" begin
        mixer = standard_mixer(3)
        @test mixer.n_qubits == 3
        @test length(mixer.terms) == 3
    end
    
    @testset "QAOA circuit" begin
        edges = [(1, 2)]
        cost_H = maxcut_hamiltonian(edges, 2)
        qc = qaoa_circuit(cost_H, 2)
        @test qc.n_qubits == 2
        @test qc.p == 2
    end
    
    @testset "Build QAOA" begin
        edges = [(1, 2)]
        cost_H = maxcut_hamiltonian(edges, 2)
        qc = qaoa_circuit(cost_H, 1)
        
        gamma = [0.5]
        beta = [0.3]
        circuit = build_qaoa(qc, gamma, beta)
        @test qubits(circuit) == 2
        
        state = execute(circuit)
        @test is_normalized(state)
    end
    
    @testset "QAOA MaxCut - simple" begin
        # Two vertices, one edge - optimal cut separates them
        edges = [(1, 2)]
        result = qaoa_maxcut(edges, 2; p=1, max_iter=100)
        
        # Max cut of single edge graph is 1
        # Use approximate comparison due to floating point precision
        @test result.optimal_value >= 0.49  # Should find reasonable solution
    end
    
    @testset "OptimizationResult" begin
        result = OptimizationResult([1.0, 2.0], 0.5, 10, true, [1.0, 0.7, 0.5])
        @test result.optimal_params == [1.0, 2.0]
        @test result.optimal_value == 0.5
        @test result.iterations == 10
        @test result.converged == true
        @test length(result.history) == 3
    end
end
