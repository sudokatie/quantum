# Quantum Phase Estimation algorithm

"""
    cphase(k::Int) -> Matrix{ComplexF64}

Controlled phase rotation gate R_k (rotation by 2pi/2^k).
"""
function cphase(k::Int)
    theta = 2 * pi / (2^k)
    ComplexF64[1 0; 0 exp(im * theta)]
end

"""
    cphase_dag(k::Int) -> Matrix{ComplexF64}

Conjugate transpose of controlled phase rotation (R_k dagger).
"""
function cphase_dag(k::Int)
    theta = 2 * pi / (2^k)
    ComplexF64[1 0; 0 exp(-im * theta)]
end

"""
    inverse_qft!(state::StateVector, qubits::Vector{Int}) -> StateVector

Apply inverse Quantum Fourier Transform in-place on specified qubits.
Qubits are 1-indexed. The first qubit in the vector is the most significant.

The inverse QFT reverses the standard QFT circuit:
- QFT: for each qubit j from 1 to n, apply H then controlled rotations, then swap
- IQFT: swap first, then for each qubit j from n to 1, apply inverse rotations then H
"""
function inverse_qft!(state::StateVector, qubits::Vector{Int})
    n = length(qubits)

    # Step 1: Swap qubits to reverse order
    for i in 1:div(n, 2)
        q1 = qubits[i]
        q2 = qubits[n - i + 1]
        state = apply(state, SWAP, [q1, q2])
    end

    # Step 2: Apply inverse QFT gates (reversed order from QFT)
    for j in n:-1:1
        target = qubits[j]

        # Apply inverse controlled phase rotations from higher-indexed qubits
        for k in n:-1:(j+1)
            control = qubits[k]
            dist = k - j + 1
            state = apply_controlled(state, cphase_dag(dist), [control], [target])
        end

        # Apply Hadamard (H is self-inverse)
        state = apply(state, H, target)
    end

    state
end

"""
    controlled_U!(state::StateVector, control::Int, target_qubits::Vector{Int},
                  U::Function, power::Int) -> StateVector

Apply U^(2^power) controlled by the control qubit.
U is a function that takes a state and target qubits and returns the transformed state.
"""
function controlled_U!(state::StateVector, control::Int, target_qubits::Vector{Int},
                       U::Function, power::Int)
    n = state.n_qubits
    dim = 2^n
    new_amps = copy(state.amplitudes)

    # Bit position for control qubit
    ctrl_pos = n - control

    # Find all basis states where control is |1>
    # For those states, we need to apply U^(2^power) to the target qubits

    # Build the subspace where control = 1
    # Extract amplitudes, apply U, put back

    # Get indices where control bit is 1
    indices_ctrl_1 = [i for i in 0:(dim-1) if (i >> ctrl_pos) & 1 == 1]

    # Create a mapping from target qubit subspace
    # We need to apply U^(2^power) to the target register
    target_dim = 2^length(target_qubits)

    # For the subspace where control=1, create a "virtual" state vector
    # over the target qubits and apply U

    # Group indices by non-target bits (excluding control which is fixed to 1)
    non_target_bits = [q for q in 1:n if q != control && !(q in target_qubits)]

    # For each configuration of non-target, non-control bits
    for other_config in 0:(2^length(non_target_bits) - 1)
        # Build the subset of amplitude indices for this configuration
        # where control=1 and non-target bits match other_config
        subset_indices = Int[]

        for target_config in 0:(target_dim - 1)
            # Reconstruct full index
            idx = 0
            # Set control bit
            idx |= (1 << ctrl_pos)

            # Set non-target bits
            for (bit_idx, q) in enumerate(non_target_bits)
                bit_val = (other_config >> (bit_idx - 1)) & 1
                pos = n - q
                idx |= (bit_val << pos)
            end

            # Set target bits
            for (bit_idx, q) in enumerate(target_qubits)
                bit_val = (target_config >> (bit_idx - 1)) & 1
                pos = n - q
                idx |= (bit_val << pos)
            end

            push!(subset_indices, idx)
        end

        # Extract amplitudes for this subspace
        sub_amps = ComplexF64[state.amplitudes[idx + 1] for idx in subset_indices]

        # Track original norm of subspace
        orig_norm = sqrt(sum(abs2, sub_amps))

        if orig_norm > 1e-15
            sub_state = StateVector(sub_amps; normalize=false)

            # Apply U^(2^power)
            for _ in 1:(2^power)
                sub_state = U(sub_state, collect(1:length(target_qubits)))
            end

            # Restore original norm (U is unitary so it preserves norm,
            # but apply() normalizes which can change the relative scale)
            current_norm = sqrt(sum(abs2, sub_state.amplitudes))
            if current_norm > 1e-15
                scale = orig_norm / current_norm
                for (i, idx) in enumerate(subset_indices)
                    new_amps[idx + 1] = sub_state.amplitudes[i] * scale
                end
            end
        end
    end

    StateVector(new_amps; normalize=false)
end

"""
    quantum_phase_estimation(eigenstate::StateVector, U::Function,
                             precision_qubits::Int) -> StateVector

Perform quantum phase estimation.

Arguments:
- eigenstate: The eigenstate |u> of the unitary U (as a StateVector)
- U: A function U(state, qubits) that applies the unitary to specified qubits
- precision_qubits: Number of qubits to use for phase precision

Returns a state where the precision register encodes the phase phi in |phi/2^n>.
The first `precision_qubits` qubits are the precision register.
The remaining qubits hold the eigenstate.

For an eigenstate |u> with U|u> = e^(2*pi*i*phi)|u>, the precision register
will approximately encode the binary representation of phi.
"""
function quantum_phase_estimation(eigenstate::StateVector, U::Function,
                                   precision_qubits::Int)
    n_eigen = num_qubits(eigenstate)
    n_total = precision_qubits + n_eigen

    if n_total > 20
        error("Total qubits ($n_total) exceeds maximum of 20")
    end

    # Create initial state: |0...0> tensor |eigenstate>
    # Precision qubits are 1:precision_qubits
    # Eigenstate qubits are (precision_qubits+1):n_total
    init_amps = zeros(ComplexF64, 2^n_total)

    # |0>^precision tensor |eigenstate>
    for (i, amp) in enumerate(eigenstate.amplitudes)
        # i-1 is the eigenstate basis index (0-indexed)
        # Full index: precision bits = 0, eigenstate bits = i-1
        full_idx = i - 1  # eigenstate bits go in lower positions
        init_amps[full_idx + 1] = amp
    end

    state = StateVector(init_amps; normalize=false)

    # Step 1: Apply Hadamard to all precision qubits
    for q in 1:precision_qubits
        state = apply(state, H, q)
    end

    # Step 2: Apply controlled-U^(2^k) operations
    # Qubit 1 controls U^(2^(n-1)), qubit 2 controls U^(2^(n-2)), etc.
    # This is the standard QPE convention
    target_qubits = collect((precision_qubits + 1):n_total)

    for k in 1:precision_qubits
        power = precision_qubits - k
        state = controlled_U!(state, k, target_qubits, U, power)
    end

    # Step 3: Apply inverse QFT to precision register
    state = inverse_qft!(state, collect(1:precision_qubits))

    state
end

"""
    estimate_phase(eigenstate::StateVector, U::Function, precision_qubits::Int;
                   shots::Int=1000) -> Float64

Estimate the phase of eigenvalue e^(2*pi*i*phi) for eigenstate of U.

Runs QPE and measures the precision register multiple times.
Returns the estimated phase as a value in [0, 1) (multiply by 2*pi for radians).
"""
function estimate_phase(eigenstate::StateVector, U::Function, precision_qubits::Int;
                        shots::Int=1000)
    # Run QPE
    state = quantum_phase_estimation(eigenstate, U, precision_qubits)

    # Sample the state
    counts = sample(state, shots)

    # Find most frequent measurement outcome
    # We only care about the precision register (first precision_qubits bits)
    n_total = num_qubits(state)
    n_eigen = n_total - precision_qubits

    # Aggregate counts by precision register value
    precision_counts = Dict{Int, Int}()
    for (outcome, count) in counts
        # Extract precision register bits (most significant bits)
        precision_val = outcome >> n_eigen
        precision_counts[precision_val] = get(precision_counts, precision_val, 0) + count
    end

    # Find the most likely precision value
    max_count = 0
    best_val = 0
    for (val, count) in precision_counts
        if count > max_count
            max_count = count
            best_val = val
        end
    end

    # Convert to phase: phi = measured_value / 2^precision_qubits
    phase = best_val / (2^precision_qubits)

    phase
end
