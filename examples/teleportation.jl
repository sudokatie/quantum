#!/usr/bin/env julia
# Quantum Teleportation Example
# Teleports an arbitrary qubit state using entanglement

using Quantum

println("=== Quantum Teleportation ===\n")

# The state to teleport: |ψ⟩ = α|0⟩ + β|1⟩
# We'll use |+⟩ = (|0⟩ + |1⟩)/√2
println("State to teleport:")
psi = apply(zero_state(1), H, 1)
show_state(psi)
println("Bloch representation:")
show_bloch(psi)

# Full teleportation circuit on 3 qubits:
# Qubit 1: Alice's qubit to teleport
# Qubit 2: Alice's half of Bell pair
# Qubit 3: Bob's half of Bell pair

println("\n--- Building teleportation circuit ---\n")

c = Circuit(3)

# Step 1: Create Bell pair between qubits 2 and 3
println("Step 1: Create Bell pair (qubits 2-3)")
add!(c, H, 2)
add!(c, CNOT, 2, 3)

# Step 2: Alice's operations
println("Step 2: Alice applies CNOT and H")
add!(c, CNOT, 1, 2)
add!(c, H, 1)

# At this point, Alice would measure qubits 1 and 2
# and send classical bits to Bob

# Step 3: Bob's corrections (controlled by Alice's measurements)
# In a real protocol, these would be conditional
# Here we show the circuit structure
println("Step 3: Conditional corrections (CNOT, CZ)")
add_controlled!(c, X, [2], [3])  # If Alice's qubit 2 measured 1
add_controlled!(c, Z, [1], [3])  # If Alice's qubit 1 measured 1

println("\nFull circuit:")
show_circuit(c)

# To properly test, we need to prepare qubit 1 in state |ψ⟩
# The circuit expects |ψ⟩|00⟩ as input

# Create initial state: |+⟩|00⟩
initial = tensor(psi, zero_state(2))
println("\nInitial state (|ψ⟩ ⊗ |00⟩):")
println("Total qubits: ", initial.n_qubits)

# Execute
final = execute(c, initial)

println("\nFinal state after teleportation:")
show_state(final)

# Note: In the full protocol, after Alice measures qubits 1,2
# and Bob applies corrections, qubit 3 will be in state |ψ⟩

println("\n--- Simplified Demonstration ---")
println("\nThe teleportation works! Bob's qubit (3) receives")
println("Alice's original state through entanglement,")
println("classical communication, and local operations.")
