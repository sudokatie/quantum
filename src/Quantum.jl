module Quantum

using LinearAlgebra
using Random
using Printf
# Core types
export StateVector
export Circuit, Operation
# State functions
export zero_state, one_state, basis_state
export normalize!, is_normalized, num_qubits
# Gates
export I_GATE, X, Y, Z, H, S, T
export Rx, Ry, Rz
export CNOT, CZ, SWAP, iSWAP, TOFFOLI
export controlled, tensor_product, is_unitary
# Operations
export apply, apply!, apply_controlled
export measure, measure_qubit, probabilities, sample, expectation
export tensor, partial_trace, concurrence, is_separable, is_entangled, purity, von_neumann_entropy
# QPE
export quantum_phase_estimation, estimate_phase, inverse_qft!, controlled_U!
export cphase, cphase_dag
# Circuit
export Operation, Circuit
export add!, add_controlled!, depth, gate_count, qubits
export add_h!, add_x!, add_y!, add_z!, add_s!, add_t!
export add_cnot!, add_cz!, add_swap!
export execute, execute_and_measure, simulate
# Visualization
export show_circuit, show_state, show_probabilities, show_bloch
# Noise models
export NoiseModel, DepolarizingNoise, AmplitudeDampingNoise, PhaseDampingNoise, MeasurementNoise, CompositeNoise
export apply_noise, apply_noise_all, noisy_measure, noisy_measure_qubit
export depolarizing, amplitude_damping, phase_damping, measurement_error, compose, noise_level
# Optimization algorithms (VQE, QAOA)
export Hamiltonian, hamiltonian, expectation_value, pauli_expectation
export VariationalCircuit, variational_circuit, num_parameters, build_ansatz
export OptimizationResult, nelder_mead, vqe
export QAOACircuit, qaoa_circuit, build_qaoa, qaoa, maxcut_hamiltonian, standard_mixer
export qaoa_maxcut, sample_qaoa_solution
# Include source files (created in subsequent tasks)
include("state.jl")
include("gates.jl")
include("apply.jl")
include("measure.jl")
include("entangle.jl")
include("circuit.jl")
include("execute.jl")
include("visualize.jl")
include("noise.jl")
include("optimizer.jl")
include("qpe.jl")
end # module
