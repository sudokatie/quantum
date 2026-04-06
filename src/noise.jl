# Noise models for quantum simulation

"""
    NoiseModel

Abstract type for noise models. Noise can be applied to quantum states
to simulate realistic quantum computers.
"""
abstract type NoiseModel end

"""
    DepolarizingNoise <: NoiseModel

Depolarizing noise channel. With probability p, replaces the qubit state
with the maximally mixed state (equivalent to randomly applying I, X, Y, or Z).

For a single qubit:
- With probability 1-p: no error
- With probability p/3: apply X
- With probability p/3: apply Y  
- With probability p/3: apply Z
"""
struct DepolarizingNoise <: NoiseModel
    p::Float64  # Error probability
    
    function DepolarizingNoise(p::Float64)
        if p < 0 || p > 1
            error("Probability must be in [0, 1], got $p")
        end
        new(p)
    end
end

"""
    AmplitudeDampingNoise <: NoiseModel

Amplitude damping channel. Models T1 relaxation (energy decay from |1> to |0>).
With probability gamma, the |1> state decays to |0>.

Kraus operators:
- K0 = [[1, 0], [0, sqrt(1-gamma)]]
- K1 = [[0, sqrt(gamma)], [0, 0]]
"""
struct AmplitudeDampingNoise <: NoiseModel
    gamma::Float64  # Decay probability
    
    function AmplitudeDampingNoise(gamma::Float64)
        if gamma < 0 || gamma > 1
            error("Gamma must be in [0, 1], got $gamma")
        end
        new(gamma)
    end
end

"""
    PhaseDampingNoise <: NoiseModel

Phase damping channel. Models T2 dephasing (loss of phase coherence).
With probability gamma, the off-diagonal elements of the density matrix decay.

Kraus operators:
- K0 = [[1, 0], [0, sqrt(1-gamma)]]
- K1 = [[0, 0], [0, sqrt(gamma)]]
"""
struct PhaseDampingNoise <: NoiseModel
    gamma::Float64  # Dephasing probability
    
    function PhaseDampingNoise(gamma::Float64)
        if gamma < 0 || gamma > 1
            error("Gamma must be in [0, 1], got $gamma")
        end
        new(gamma)
    end
end

"""
    MeasurementNoise <: NoiseModel

Measurement error model. With probability p, the measurement outcome is flipped.
- |0> measured as |1> with probability p
- |1> measured as |0> with probability p
"""
struct MeasurementNoise <: NoiseModel
    p::Float64  # Flip probability
    
    function MeasurementNoise(p::Float64)
        if p < 0 || p > 0.5
            error("Flip probability must be in [0, 0.5], got $p")
        end
        new(p)
    end
end

"""
    CompositeNoise <: NoiseModel

Combines multiple noise models.
"""
struct CompositeNoise <: NoiseModel
    models::Vector{NoiseModel}
end

"""
    apply_noise(state::StateVector, noise::NoiseModel, qubit::Int) -> StateVector

Apply noise to a specific qubit in the state.
"""
function apply_noise(state::StateVector, noise::DepolarizingNoise, qubit::Int)
    if qubit < 1 || qubit > state.n_qubits
        error("Qubit index out of range: $qubit (state has $(state.n_qubits) qubits)")
    end
    
    if noise.p == 0
        return state
    end
    
    # Randomly select error
    r = rand()
    if r < 1 - noise.p
        # No error
        return state
    elseif r < 1 - noise.p + noise.p/3
        # Apply X
        return apply(state, X, qubit)
    elseif r < 1 - noise.p + 2*noise.p/3
        # Apply Y
        return apply(state, Y, qubit)
    else
        # Apply Z
        return apply(state, Z, qubit)
    end
end

function apply_noise(state::StateVector, noise::AmplitudeDampingNoise, qubit::Int)
    if qubit < 1 || qubit > state.n_qubits
        error("Qubit index out of range: $qubit (state has $(state.n_qubits) qubits)")
    end
    
    if noise.gamma == 0
        return state
    end
    
    # Kraus operators for amplitude damping
    # K0 = [[1, 0], [0, sqrt(1-gamma)]]
    # K1 = [[0, sqrt(gamma)], [0, 0]]
    
    n = state.n_qubits
    dim = 2^n
    new_amps = zeros(ComplexF64, dim)
    
    sqrt_1_minus_gamma = sqrt(1 - noise.gamma)
    sqrt_gamma = sqrt(noise.gamma)
    
    # Apply Kraus operators
    for i in 0:(dim-1)
        # Check if qubit is 0 or 1 in this basis state
        qubit_val = (i >> (n - qubit)) & 1
        
        if qubit_val == 0
            # |0> contribution from K0|0> = |0>
            new_amps[i+1] += state.amplitudes[i+1]
        else
            # |1> contribution
            # K0|1> = sqrt(1-gamma)|1>
            new_amps[i+1] += sqrt_1_minus_gamma * state.amplitudes[i+1]
            
            # K1|1> = sqrt(gamma)|0>
            # Flip the qubit bit to get the |0> state
            j = i ⊻ (1 << (n - qubit))
            new_amps[j+1] += sqrt_gamma * state.amplitudes[i+1]
        end
    end
    
    StateVector(new_amps)
end

function apply_noise(state::StateVector, noise::PhaseDampingNoise, qubit::Int)
    if qubit < 1 || qubit > state.n_qubits
        error("Qubit index out of range: $qubit (state has $(state.n_qubits) qubits)")
    end
    
    if noise.gamma == 0
        return state
    end
    
    # Kraus operators for phase damping
    # K0 = [[1, 0], [0, sqrt(1-gamma)]]
    # K1 = [[0, 0], [0, sqrt(gamma)]]
    
    n = state.n_qubits
    dim = 2^n
    new_amps = zeros(ComplexF64, dim)
    
    sqrt_1_minus_gamma = sqrt(1 - noise.gamma)
    
    # Apply K0 (dominates for small gamma)
    for i in 0:(dim-1)
        qubit_val = (i >> (n - qubit)) & 1
        
        if qubit_val == 0
            new_amps[i+1] += state.amplitudes[i+1]
        else
            new_amps[i+1] += sqrt_1_minus_gamma * state.amplitudes[i+1]
        end
    end
    
    # K1 only affects populations, not coherences for pure states
    # For a full density matrix simulation, we'd need to track the density matrix
    # For state vector simulation, phase damping is approximated by the above
    
    StateVector(new_amps)
end

function apply_noise(state::StateVector, noise::CompositeNoise, qubit::Int)
    result = state
    for model in noise.models
        result = apply_noise(result, model, qubit)
    end
    return result
end

"""
    apply_noise_all(state::StateVector, noise::NoiseModel) -> StateVector

Apply noise to all qubits in the state.
"""
function apply_noise_all(state::StateVector, noise::NoiseModel)
    result = state
    for q in 1:state.n_qubits
        result = apply_noise(result, noise, q)
    end
    return result
end

"""
    noisy_measure(state::StateVector, noise::MeasurementNoise) -> (Int, StateVector)

Measure all qubits with measurement noise. Returns the (possibly noisy) outcome
and the collapsed state.
"""
function noisy_measure(state::StateVector, noise::MeasurementNoise)
    # First, do ideal measurement
    outcome, collapsed = measure(state)
    
    if noise.p == 0
        return (outcome, collapsed)
    end
    
    # Convert outcome to bits, apply noise, convert back
    n = state.n_qubits
    bits = [(outcome >> (n - i)) & 1 for i in 1:n]
    
    # Apply measurement errors
    for i in 1:n
        if rand() < noise.p
            bits[i] = 1 - bits[i]  # Flip 0->1 or 1->0
        end
    end
    
    # Convert back to integer
    noisy_outcome = sum(bits[i] << (n - i) for i in 1:n)
    
    return (noisy_outcome, collapsed)
end

"""
    noisy_measure_qubit(state::StateVector, qubit::Int, noise::MeasurementNoise) -> (Int, StateVector)

Measure a single qubit with measurement noise.
"""
function noisy_measure_qubit(state::StateVector, qubit::Int, noise::MeasurementNoise)
    # First, do ideal measurement
    outcome, new_state = measure_qubit(state, qubit)
    
    if noise.p > 0 && rand() < noise.p
        outcome = 1 - outcome  # Flip outcome
    end
    
    return outcome, new_state
end

# Convenience constructors

"""
    depolarizing(p::Float64) -> DepolarizingNoise

Create a depolarizing noise channel with error probability p.
"""
depolarizing(p::Float64) = DepolarizingNoise(p)

"""
    amplitude_damping(gamma::Float64) -> AmplitudeDampingNoise

Create an amplitude damping channel with decay probability gamma.
"""
amplitude_damping(gamma::Float64) = AmplitudeDampingNoise(gamma)

"""
    phase_damping(gamma::Float64) -> PhaseDampingNoise

Create a phase damping channel with dephasing probability gamma.
"""
phase_damping(gamma::Float64) = PhaseDampingNoise(gamma)

"""
    measurement_error(p::Float64) -> MeasurementNoise

Create a measurement error model with flip probability p.
"""
measurement_error(p::Float64) = MeasurementNoise(p)

"""
    compose(models::NoiseModel...) -> CompositeNoise

Combine multiple noise models.
"""
compose(models::NoiseModel...) = CompositeNoise(collect(models))

# Noise level presets

"""
    noise_level(level::Symbol) -> NoiseModel

Get a preset noise model.
- :none - No noise
- :low - Low noise (p=0.001)
- :medium - Medium noise (p=0.01)
- :high - High noise (p=0.1)
- :ibm_like - Realistic IBM-like noise
"""
function noise_level(level::Symbol)
    if level == :none
        return DepolarizingNoise(0.0)
    elseif level == :low
        return DepolarizingNoise(0.001)
    elseif level == :medium
        return DepolarizingNoise(0.01)
    elseif level == :high
        return DepolarizingNoise(0.1)
    elseif level == :ibm_like
        # Typical IBM quantum computer noise levels
        return CompositeNoise([
            DepolarizingNoise(0.001),      # Single-qubit gate error
            AmplitudeDampingNoise(0.0001), # T1 decay
            PhaseDampingNoise(0.0002)      # T2 dephasing
        ])
    else
        error("Unknown noise level: $level")
    end
end
