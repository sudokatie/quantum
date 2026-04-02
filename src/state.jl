# State vector representation

"""
    StateVector

Represents a quantum state as a vector of complex amplitudes.
For n qubits, has 2^n amplitudes in computational basis order.
"""
struct StateVector
    amplitudes::Vector{ComplexF64}
    n_qubits::Int
    
    function StateVector(amplitudes::Vector{ComplexF64}; normalize::Bool=true)
        n = length(amplitudes)
        # Check power of 2
        n_qubits = round(Int, log2(n))
        if 2^n_qubits != n
            error("Amplitude vector length must be power of 2, got $n")
        end
        # Enforce maximum 20 qubits per spec
        if n_qubits > 20
            error("Maximum 20 qubits supported, got $n_qubits")
        end
        # Auto-normalize per spec
        if normalize
            norm_factor = sqrt(sum(abs2, amplitudes))
            if norm_factor > 1e-15
                normalized = amplitudes ./ norm_factor
                new(normalized, n_qubits)
            else
                new(amplitudes, n_qubits)
            end
        else
            new(amplitudes, n_qubits)
        end
    end
end

"""
    zero_state(n::Int) -> StateVector

Create the |0>^n state (all qubits in |0>).
"""
function zero_state(n::Int)
    if n > 20
        error("Maximum 20 qubits supported, got $n")
    end
    amps = zeros(ComplexF64, 2^n)
    amps[1] = 1.0
    StateVector(amps; normalize=false)
end

"""
    one_state(n::Int) -> StateVector

Create the |1>^n state (all qubits in |1>).
"""
function one_state(n::Int)
    if n > 20
        error("Maximum 20 qubits supported, got $n")
    end
    amps = zeros(ComplexF64, 2^n)
    amps[end] = 1.0
    StateVector(amps; normalize=false)
end

"""
    basis_state(n::Int, index::Int) -> StateVector

Create computational basis state |index> for n qubits.
Index is 0-based (0 to 2^n - 1).
"""
function basis_state(n::Int, index::Int)
    if n > 20
        error("Maximum 20 qubits supported, got $n")
    end
    if index < 0 || index >= 2^n
        error("Index $index out of range for $n qubits")
    end
    amps = zeros(ComplexF64, 2^n)
    amps[index + 1] = 1.0
    StateVector(amps; normalize=false)
end

"""
    num_qubits(state::StateVector) -> Int

Return the number of qubits in the state.
"""
num_qubits(state::StateVector) = state.n_qubits

"""
    is_normalized(state::StateVector; tol=1e-10) -> Bool

Check if state vector is normalized (|amplitudes|^2 sum to 1).
"""
function is_normalized(state::StateVector; tol=1e-10)
    total = sum(abs2, state.amplitudes)
    abs(total - 1.0) < tol
end

"""
    normalize!(state::StateVector) -> StateVector

Normalize state vector in place.
"""
function normalize!(state::StateVector)
    norm = sqrt(sum(abs2, state.amplitudes))
    state.amplitudes ./= norm
    state
end

# Equality
Base.:(==)(a::StateVector, b::StateVector) = a.amplitudes == b.amplitudes

# Copy
Base.copy(s::StateVector) = StateVector(copy(s.amplitudes))

# Length
Base.length(s::StateVector) = length(s.amplitudes)

# Indexing
Base.getindex(s::StateVector, i) = s.amplitudes[i]
