@testset "Quantum Phase Estimation" begin
    @testset "cphase gates" begin
        # R_1 should be Z gate (rotation by pi)
        r1 = cphase(1)
        @test r1[1, 1] == 1.0
        @test abs(r1[2, 2] - (-1.0)) < 1e-10  # e^(i*pi) = -1

        # R_2 should be S gate (rotation by pi/2)
        r2 = cphase(2)
        @test r2[1, 1] == 1.0
        @test abs(r2[2, 2] - im) < 1e-10  # e^(i*pi/2) = i

        # Conjugate should be inverse
        r2_dag = cphase_dag(2)
        @test abs(r2_dag[2, 2] - conj(r2[2, 2])) < 1e-10
    end

    @testset "inverse QFT" begin
        # IQFT of |00> gives uniform superposition (same as QFT)
        state = zero_state(2)
        result = inverse_qft!(state, [1, 2])
        # Should be close to uniform: each amplitude ~ 0.5
        @test abs(abs(result.amplitudes[1]) - 0.5) < 0.1
        @test abs(abs(result.amplitudes[2]) - 0.5) < 0.1

        # IQFT of uniform superposition gives |00>
        # This is the key property: IQFT(QFT(|00>)) = |00>
        state = zero_state(2)
        state = apply(state, H, 1)
        state = apply(state, H, 2)
        result = inverse_qft!(state, [1, 2])
        @test abs(result.amplitudes[1]) > 0.99
    end

    @testset "controlled_U!" begin
        # Test with a simple phase gate as U
        # U = Z gate, eigenstate |1> has eigenvalue -1 = e^(i*pi)
        function apply_z(state, qubits)
            apply(state, Z, qubits[1])
        end

        # Create |+> tensor |1> (control in superposition, target in eigenstate)
        state = zero_state(2)
        state = apply(state, H, 1)       # Control in |+>
        state = apply(state, X, 2)       # Target in |1>

        # Apply controlled-Z (Z^1)
        result = controlled_U!(state, 1, [2], apply_z, 0)

        # When control is |1>, Z acts on |1> giving -|1>
        # |+>|1> = (|0> + |1>)|1>/sqrt(2) = (|01> + |11>)/sqrt(2)
        # After controlled-Z: (|01> - |11>)/sqrt(2) = |->|1>
        # |01> has index 1, |11> has index 3 (Julia 1-indexed: 2 and 4)
        @test abs(result.amplitudes[2]) > 0.7    # |01> amplitude
        @test abs(result.amplitudes[4]) > 0.7    # |11> amplitude (magnitude same)
        # Check opposite phases (one positive, one negative)
        @test real(result.amplitudes[2]) * real(result.amplitudes[4]) < 0
    end

    @testset "QPE with T gate" begin
        # T gate has eigenvalue e^(i*pi/4) for eigenstate |1>
        # Phase = 1/8 (since e^(2*pi*i*1/8) = e^(i*pi/4))
        function apply_t(state, qubits)
            apply(state, T, qubits[1])
        end

        # Eigenstate |1>
        eigenstate = one_state(1)

        # Run QPE with 3 precision qubits
        result = quantum_phase_estimation(eigenstate, apply_t, 3)

        # The precision register should measure 1 (binary 001 = 1/8)
        # State should have high amplitude at |001>|1> = |0011> = index 3 (0-indexed)
        # With 4 total qubits: precision (3) + eigenstate (1)
        # |001>|1> in big-endian: qubit1=0, qubit2=0, qubit3=1, qubit4=1
        # That's binary 0011 = 3

        @test num_qubits(result) == 4
        @test abs(result.amplitudes[4])^2 > 0.5  # Should have high probability
    end

    @testset "estimate_phase convergence" begin
        # T gate: eigenvalue e^(i*pi/4), phase = 1/8
        function apply_t(state, qubits)
            apply(state, T, qubits[1])
        end

        eigenstate = one_state(1)

        # Test with increasing precision
        phase_3 = estimate_phase(eigenstate, apply_t, 3; shots=1000)
        phase_4 = estimate_phase(eigenstate, apply_t, 4; shots=1000)
        phase_5 = estimate_phase(eigenstate, apply_t, 5; shots=1000)

        expected = 1/8  # 0.125

        # With 3 qubits, should get exactly 1/8 (since 1/8 = 1/2^3)
        @test abs(phase_3 - expected) < 0.02

        # Higher precision should be at least as good
        @test abs(phase_4 - expected) < 0.1
        @test abs(phase_5 - expected) < 0.1
    end

    @testset "estimate_phase with S gate" begin
        # S gate has eigenvalue i = e^(i*pi/2) for |1>
        # Phase = 1/4 (since e^(2*pi*i*1/4) = e^(i*pi/2) = i)
        function apply_s(state, qubits)
            apply(state, S, qubits[1])
        end

        eigenstate = one_state(1)

        # Use more precision qubits for better accuracy
        phase = estimate_phase(eigenstate, apply_s, 4; shots=1000)
        expected = 1/4  # 0.25

        @test abs(phase - expected) < 0.05
    end

    @testset "estimate_phase with Z gate" begin
        # Z gate has eigenvalue -1 = e^(i*pi) for |1>
        # Phase = 1/2
        function apply_z(state, qubits)
            apply(state, Z, qubits[1])
        end

        eigenstate = one_state(1)

        phase = estimate_phase(eigenstate, apply_z, 3; shots=1000)
        expected = 1/2  # 0.5

        @test abs(phase - expected) < 0.01
    end

    @testset "QPE with 2-qubit eigenstate" begin
        # Test with a 2-qubit eigenstate
        # Use controlled-Z which has eigenvalue -1 for |11>
        function apply_cz(state, qubits)
            apply(state, CZ, qubits)
        end

        # Eigenstate |11>
        eigenstate = basis_state(2, 3)  # |11>

        # CZ|11> = -|11>, so phase = 1/2
        phase = estimate_phase(eigenstate, apply_cz, 3; shots=1000)

        @test abs(phase - 0.5) < 0.05
    end
end
