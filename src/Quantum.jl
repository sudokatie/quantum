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
export CNOT, CZ, SWAP, iSWAP
export controlled, tensor_product, is_unitary

# Operations
export apply, apply!
export measure, measure_qubit, probabilities, sample, expectation
export tensor, partial_trace, concurrence, is_separable, purity, von_neumann_entropy

# Circuit
export Operation, Circuit
export add!, add_controlled!, depth, gate_count, qubits
export add_h!, add_x!, add_y!, add_z!, add_s!, add_t!
export add_cnot!, add_cz!, add_swap!
export execute, execute_and_measure, simulate

# Visualization
export show_circuit, show_state, show_probabilities, show_bloch

# Include source files (created in subsequent tasks)
include("state.jl")
include("gates.jl")
include("apply.jl")
include("measure.jl")
include("entangle.jl")
include("circuit.jl")
include("execute.jl")
include("visualize.jl")

end # module
