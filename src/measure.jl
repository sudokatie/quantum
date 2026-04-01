# Measurement operations

"""
    probabilities(state::StateVector) -> Vector{Float64}

Get the probability distribution for measuring each computational basis state.
"""
function probabilities(state::StateVector)
    [abs2(a) for a in state.amplitudes]
end

"""
    measure(state::StateVector) -> Tuple{Int, StateVector}

Measure all qubits, collapsing the state.
Returns (result, collapsed_state) where result is the measured basis state (0-indexed).
"""
function measure(state::StateVector)
    probs = probabilities(state)
    
    # Sample from the distribution
    r = rand()
    cumulative = 0.0
    result = 0
    
    for (i, p) in enumerate(probs)
        cumulative += p
        if r < cumulative
            result = i - 1  # Convert to 0-indexed
            break
        end
    end
    
    # Collapse to the measured state
    collapsed = basis_state(state.n_qubits, result)
    
    (result, collapsed)
end

"""
    measure_qubit(state::StateVector, qubit::Int) -> Tuple{Int, StateVector}

Measure a single qubit, partially collapsing the state.
Returns (bit, collapsed_state) where bit is 0 or 1.
"""
function measure_qubit(state::StateVector, qubit::Int)
    @assert 1 <= qubit <= state.n_qubits "Qubit index out of range"
    
    n = state.n_qubits
    t = n - qubit  # Bit position (qubit 1 = MSB)
    
    # Calculate probability of measuring 0
    prob_0 = 0.0
    for i in 0:(2^n - 1)
        if (i >> t) & 1 == 0
            prob_0 += abs2(state.amplitudes[i + 1])
        end
    end
    
    # Sample the result
    result = rand() < prob_0 ? 0 : 1
    
    # Collapse the state
    new_amps = zeros(ComplexF64, 2^n)
    norm_factor = 0.0
    
    for i in 0:(2^n - 1)
        bit_val = (i >> t) & 1
        if bit_val == result
            new_amps[i + 1] = state.amplitudes[i + 1]
            norm_factor += abs2(state.amplitudes[i + 1])
        end
    end
    
    # Normalize
    norm_factor = sqrt(norm_factor)
    new_amps ./= norm_factor
    
    (result, StateVector(new_amps))
end

"""
    sample(state::StateVector, shots::Int) -> Dict{Int, Int}

Sample the state multiple times without modifying it.
Returns a dictionary mapping basis state (0-indexed) to count.
"""
function sample(state::StateVector, shots::Int)
    probs = probabilities(state)
    counts = Dict{Int, Int}()
    
    for _ in 1:shots
        r = rand()
        cumulative = 0.0
        
        for (i, p) in enumerate(probs)
            cumulative += p
            if r < cumulative
                result = i - 1
                counts[result] = get(counts, result, 0) + 1
                break
            end
        end
    end
    
    counts
end

"""
    expectation(state::StateVector, observable::Matrix{ComplexF64}) -> Float64

Calculate the expectation value of an observable.
observable must be a Hermitian matrix of size 2^n x 2^n.
"""
function expectation(state::StateVector, observable::Matrix{ComplexF64})
    @assert size(observable) == (2^state.n_qubits, 2^state.n_qubits) "Observable size mismatch"
    
    # <psi|O|psi> = conj(amplitudes) * O * amplitudes
    result = state.amplitudes' * observable * state.amplitudes
    real(result)  # Should be real for Hermitian observable
end
