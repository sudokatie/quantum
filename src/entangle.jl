# Entanglement operations

"""
    tensor(a::StateVector, b::StateVector) -> StateVector

Compute tensor product of two states: |a> ⊗ |b>.
Result has a.n_qubits + b.n_qubits qubits.
"""
function tensor(a::StateVector, b::StateVector)
    # Kronecker product of amplitudes
    new_amps = kron(a.amplitudes, b.amplitudes)
    StateVector(new_amps)
end

"""
    partial_trace(state::StateVector, keep_qubits::Vector{Int}) -> Matrix{ComplexF64}

Trace out all qubits except those in keep_qubits.
Returns the reduced density matrix for the kept subsystem.
"""
function partial_trace(state::StateVector, keep_qubits::Vector{Int})
    n = state.n_qubits
    n_keep = length(keep_qubits)
    
    @assert all(1 .<= keep_qubits .<= n) "Qubit indices out of range"
    @assert length(unique(keep_qubits)) == n_keep "Duplicate qubit indices"
    
    # Dimensions
    dim_keep = 2^n_keep
    dim_trace = 2^(n - n_keep)
    
    # Result density matrix
    rho = zeros(ComplexF64, dim_keep, dim_keep)
    
    # Convert keep_qubits to bit positions (qubit 1 = MSB)
    keep_bits = sort(n .- keep_qubits, rev=true)
    trace_bits = setdiff(0:(n-1), n .- keep_qubits)
    
    for i in 0:(dim_keep - 1)
        for j in 0:(dim_keep - 1)
            # Sum over traced-out indices
            for k in 0:(dim_trace - 1)
                # Construct full indices
                full_i = reassemble_index(i, k, keep_bits, trace_bits, n)
                full_j = reassemble_index(j, k, keep_bits, trace_bits, n)
                
                rho[i + 1, j + 1] += state.amplitudes[full_i + 1] * 
                                     conj(state.amplitudes[full_j + 1])
            end
        end
    end
    
    rho
end

"""
Reassemble a full index from kept and traced indices.
"""
function reassemble_index(keep_idx::Int, trace_idx::Int, 
                          keep_bits::Vector{Int}, trace_bits::Vector{Int}, n::Int)
    result = 0
    
    # Fill in kept bits
    for (i, pos) in enumerate(keep_bits)
        bit = (keep_idx >> (length(keep_bits) - i)) & 1
        result |= (bit << pos)
    end
    
    # Fill in traced bits
    for (i, pos) in enumerate(trace_bits)
        bit = (trace_idx >> (length(trace_bits) - i)) & 1
        result |= (bit << pos)
    end
    
    result
end

"""
    purity(rho::Matrix{ComplexF64}) -> Float64

Calculate purity of a density matrix: Tr(ρ²).
Purity = 1 for pure states, < 1 for mixed states.
"""
function purity(rho::Matrix{ComplexF64})
    real(tr(rho * rho))
end

"""
    concurrence(state::StateVector) -> Float64

Calculate concurrence for a 2-qubit state.
Concurrence = 0 for separable states, 1 for maximally entangled.
"""
function concurrence(state::StateVector)
    @assert state.n_qubits == 2 "Concurrence only defined for 2-qubit states"
    
    # Build density matrix
    rho = state.amplitudes * state.amplitudes'
    
    # Pauli Y tensor Y
    sigma_y = ComplexF64[0 -im; im 0]
    yy = kron(sigma_y, sigma_y)
    
    # R = rho * (Y⊗Y) * conj(rho) * (Y⊗Y)
    rho_tilde = yy * conj(rho) * yy
    R = rho * rho_tilde
    
    # Eigenvalues of R
    eigvals_R = eigvals(R)
    
    # Take square roots of eigenvalues (may need to handle numerical issues)
    sqrt_eigvals = sort([sqrt(max(0, real(e))) for e in eigvals_R], rev=true)
    
    # Concurrence = max(0, λ1 - λ2 - λ3 - λ4)
    C = sqrt_eigvals[1] - sqrt_eigvals[2] - sqrt_eigvals[3] - sqrt_eigvals[4]
    max(0.0, C)
end

"""
    is_separable(state::StateVector; tol=1e-10) -> Bool

Check if a 2-qubit state is separable (not entangled).
Uses concurrence: separable iff concurrence ≈ 0.
"""
function is_separable(state::StateVector; tol=1e-6)
    @assert state.n_qubits == 2 "is_separable only implemented for 2-qubit states"
    concurrence(state) < tol
end

"""
    von_neumann_entropy(rho::Matrix{ComplexF64}) -> Float64

Calculate von Neumann entropy: S = -Tr(ρ log₂ ρ).
"""
function von_neumann_entropy(rho::Matrix{ComplexF64})
    eigvals_rho = eigvals(rho)
    entropy = 0.0
    
    for e in eigvals_rho
        p = real(e)
        if p > 1e-15
            entropy -= p * log2(p)
        end
    end
    
    entropy
end
