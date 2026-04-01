# Quantum

A quantum computing simulator in Julia. Implements qubit state vectors, quantum gates, measurement, entanglement, and circuit visualization.

## Features

- **State Vectors**: Represent quantum states as complex amplitude vectors
- **Quantum Gates**: Pauli (X, Y, Z), Hadamard, Phase (S, T), CNOT, SWAP, and rotations
- **Measurement**: Born rule probabilities, partial measurement, repeated sampling
- **Entanglement**: Tensor products, partial trace, concurrence, entropy
- **Circuits**: Compose gates into circuits, execute on states
- **Visualization**: ASCII circuit diagrams and state display

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/quantum")
```

## Quick Start

```julia
using Quantum

# Create a Bell state |Φ+⟩ = (|00⟩ + |11⟩) / √2
c = Circuit(2)
add!(c, H, 1)
add!(c, CNOT, 1, 2)

# View the circuit
show_circuit(c)
# q1: [H]─●─
# q2: ───⊕─

# Execute and show results
state = execute(c)
show_state(state)
# |00⟩: 0.7071
# |11⟩: 0.7071

# Sample measurements
counts = simulate(c, 1000)
# Dict(0 => 512, 3 => 488)
```

## Using Circuits

```julia
# Build a circuit
c = Circuit(3)
add!(c, H, 1)           # Hadamard on qubit 1
add!(c, CNOT, 1, 2)     # CNOT: control 1, target 2
add!(c, X, 3)           # Pauli-X on qubit 3

# Convenience methods
add_h!(c, 1)
add_x!(c, 2)
add_cnot!(c, 1, 3)

# Circuit info
gate_count(c)  # Number of operations
depth(c)       # Circuit depth
```

## Available Gates

### Single-Qubit
- `I_GATE` - Identity
- `X`, `Y`, `Z` - Pauli gates
- `H` - Hadamard
- `S`, `T` - Phase gates
- `Rx(θ)`, `Ry(θ)`, `Rz(θ)` - Rotation gates

### Two-Qubit
- `CNOT` - Controlled-NOT
- `CZ` - Controlled-Z
- `SWAP` - Swap qubits
- `iSWAP` - iSWAP gate

### Controlled Gates
```julia
# Create controlled version of any single-qubit gate
add_controlled!(circuit, X, [1], [2])  # Controlled-X
add_controlled!(circuit, Z, [1, 2], [3])  # Toffoli-like
```

## Measurement

```julia
state = execute(circuit)

# Get probability distribution
probs = probabilities(state)

# Measure all qubits (collapses state)
result, collapsed = measure(state)

# Measure single qubit
bit, collapsed = measure_qubit(state, 1)

# Sample without collapse
counts = sample(state, 1000)

# Expectation value
expectation(state, Z)
```

## Entanglement

```julia
# Tensor product of states
combined = tensor(state_a, state_b)

# Check if entangled
is_separable(state)  # false for entangled states

# Entanglement measure (2 qubits)
concurrence(state)  # 0 = separable, 1 = maximally entangled

# Partial trace (reduced density matrix)
rho = partial_trace(state, [1])  # Keep qubit 1

# Von Neumann entropy
von_neumann_entropy(rho)
```

## Visualization

```julia
# Circuit diagram
show_circuit(circuit)

# State in ket notation
show_state(state)

# Probability histogram
show_probabilities(state)

# Bloch sphere (single qubit)
show_bloch(state)
```

## Examples

See `examples/` for:
- `bell_state.jl` - Bell state preparation and measurement
- `teleportation.jl` - Quantum teleportation protocol
- `grover.jl` - Grover's search algorithm
- `qft.jl` - Quantum Fourier Transform

Run an example:
```bash
julia --project=. examples/bell_state.jl
```

## API Reference

### States
- `zero_state(n)` - Create |0⟩^⊗n
- `one_state(n)` - Create |1⟩^⊗n
- `basis_state(n, i)` - Create |i⟩
- `is_normalized(state)` - Check normalization

### Operations
- `apply(state, gate, target)` - Apply gate
- `apply(state, gate, [targets])` - Multi-qubit gate

### Circuits
- `Circuit(n)` - Create n-qubit circuit
- `add!(circuit, gate, targets...)` - Add operation
- `execute(circuit)` - Run from |0⟩^n
- `simulate(circuit, shots)` - Sample results

## Limitations

- Maximum ~20 qubits (2^20 amplitudes ≈ 1M complex numbers)
- Ideal gates only (no noise model)
- No circuit optimization passes

## Running Tests

```bash
julia --project=. test/runtests.jl
```

## License

MIT
