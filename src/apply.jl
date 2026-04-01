# Gate application to state vectors

"""
    apply(state::StateVector, gate::Matrix{ComplexF64}, target::Int) -> StateVector

Apply a single-qubit gate to the specified qubit (1-indexed).
Uses efficient algorithm that avoids materializing full matrix.
"""
function apply(state::StateVector, gate::Matrix{ComplexF64}, target::Int)
    @assert size(gate) == (2, 2) "Gate must be 2x2 for single qubit"
    @assert 1 <= target <= state.n_qubits "Target qubit out of range"
    
    new_amps = copy(state.amplitudes)
    n = state.n_qubits
    
    # Iterate through basis states, processing pairs
    # target is 1-indexed, convert to bit position
    # Qubit 1 is MSB (leftmost), so bit position = n - target
    t = n - target
    
    for i in 0:(2^n - 1)
        # Only process when target bit is 0 (to avoid double processing)
        if (i >> t) & 1 == 0
            # j is the partner state with target bit flipped
            j = i | (1 << t)
            
            # Get current amplitudes (1-indexed)
            a0 = state.amplitudes[i + 1]
            a1 = state.amplitudes[j + 1]
            
            # Apply 2x2 gate
            new_amps[i + 1] = gate[1,1] * a0 + gate[1,2] * a1
            new_amps[j + 1] = gate[2,1] * a0 + gate[2,2] * a1
        end
    end
    
    StateVector(new_amps)
end

"""
    apply(state::StateVector, gate::Matrix{ComplexF64}, targets::Vector{Int}) -> StateVector

Apply a multi-qubit gate to specified qubits.
For 2-qubit gates, targets[1] is control (if applicable), targets[2] is target.
"""
function apply(state::StateVector, gate::Matrix{ComplexF64}, targets::Vector{Int})
    n_gate_qubits = round(Int, log2(size(gate, 1)))
    @assert length(targets) == n_gate_qubits "Number of targets must match gate size"
    @assert all(1 .<= targets .<= state.n_qubits) "Target qubits out of range"
    @assert length(unique(targets)) == length(targets) "Duplicate target qubits"
    
    if n_gate_qubits == 1
        return apply(state, gate, targets[1])
    end
    
    # For multi-qubit gates, use full matrix multiplication approach
    # Build the full operator matrix
    full_gate = build_full_gate(gate, targets, state.n_qubits)
    new_amps = full_gate * state.amplitudes
    StateVector(new_amps)
end

"""
    build_full_gate(gate::Matrix, targets::Vector{Int}, n_qubits::Int) -> Matrix

Build the full 2^n x 2^n matrix for a gate acting on specific qubits.
"""
function build_full_gate(gate::Matrix{ComplexF64}, targets::Vector{Int}, n_qubits::Int)
    dim = 2^n_qubits
    result = zeros(ComplexF64, dim, dim)
    
    gate_dim = size(gate, 1)
    n_gate_qubits = round(Int, log2(gate_dim))
    
    # Convert to bit positions (qubit 1 = MSB)
    t = n_qubits .- targets
    
    for i in 0:(dim - 1)
        for j in 0:(dim - 1)
            # Extract the bits at target positions for both indices
            i_gate_bits = extract_bits(i, t)
            j_gate_bits = extract_bits(j, t)
            
            # Check if non-target bits match
            i_other = clear_bits(i, t)
            j_other = clear_bits(j, t)
            
            if i_other == j_other
                # Get the matrix element from gate
                result[i + 1, j + 1] = gate[i_gate_bits + 1, j_gate_bits + 1]
            end
        end
    end
    
    result
end

"""
Extract bits at specified positions and combine into single integer.
targets[1] maps to MSB of result, targets[end] maps to LSB.
"""
function extract_bits(n::Int, positions::Vector{Int})
    result = 0
    n_bits = length(positions)
    for (idx, pos) in enumerate(positions)
        bit = (n >> pos) & 1
        # Map to reversed position: idx=1 -> MSB, idx=n -> LSB
        result |= (bit << (n_bits - idx))
    end
    result
end

"""
Clear bits at specified positions.
"""
function clear_bits(n::Int, positions::Vector{Int})
    result = n
    for pos in positions
        result &= ~(1 << pos)
    end
    result
end

"""
    apply!(state::StateVector, gate::Matrix{ComplexF64}, target::Int) -> StateVector

Apply gate in-place (modifies state).
"""
function apply!(state::StateVector, gate::Matrix{ComplexF64}, target::Int)
    @assert size(gate) == (2, 2) "Gate must be 2x2 for single qubit"
    @assert 1 <= target <= state.n_qubits "Target qubit out of range"
    
    n = state.n_qubits
    t = n - target
    
    for i in 0:(2^n - 1)
        if (i >> t) & 1 == 0
            j = i | (1 << t)
            
            a0 = state.amplitudes[i + 1]
            a1 = state.amplitudes[j + 1]
            
            state.amplitudes[i + 1] = gate[1,1] * a0 + gate[1,2] * a1
            state.amplitudes[j + 1] = gate[2,1] * a0 + gate[2,2] * a1
        end
    end
    
    state
end
