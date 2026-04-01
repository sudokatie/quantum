#!/usr/bin/env julia
# Grover's Search Algorithm
# Searches for a marked item in an unstructured database

using Quantum

println("=== Grover's Search Algorithm ===\n")

# For N=4 items, we need n=2 qubits
# We'll search for item |11⟩ (index 3)
N = 4
n = 2
target = 3

println("Searching for |$(lpad(string(target, base=2), n, '0'))⟩ among $N items")
println("Optimal iterations: ", round(Int, pi/4 * sqrt(N)))

# Build Grover circuit
c = Circuit(n)

# Step 1: Initialize superposition
println("\nStep 1: Create uniform superposition")
for i in 1:n
    add!(c, H, i)
end

# Step 2: Grover iteration (repeat sqrt(N)/4 ≈ 1 time for N=4)
println("Step 2: Apply Grover iteration")

# Oracle: flip sign of target state |11⟩
# This is a CZ gate for target=|11⟩
add!(c, CZ, 1, 2)

# Diffusion operator: 2|s⟩⟨s| - I
# Implemented as: H⊗H · (2|0⟩⟨0| - I) · H⊗H
for i in 1:n
    add!(c, H, i)
end

# 2|0⟩⟨0| - I = -Z⊗Z (with signs)
for i in 1:n
    add!(c, X, i)
end
add!(c, CZ, 1, 2)
for i in 1:n
    add!(c, X, i)
end

for i in 1:n
    add!(c, H, i)
end

println("\nCircuit:")
show_circuit(c)

println("\nGate count: ", gate_count(c))
println("Circuit depth: ", depth(c))

# Execute
state = execute(c)

println("\n--- Results ---\n")
println("State after Grover iteration:")
show_state(state)

println("\nProbabilities:")
show_probabilities(state)

# Measure
println("\nSampling 1000 measurements:")
counts = simulate(c, 1000)
for (k, v) in sort(collect(counts))
    basis = lpad(string(k, base=2), n, '0')
    marker = k == target ? " ← TARGET" : ""
    println("  |$basis⟩: $v times$marker")
end

# Calculate success probability
target_prob = probabilities(state)[target + 1]
println("\nSuccess probability: $(round(target_prob * 100, digits=1))%")
println("Classical probability: $(round(100/N, digits=1))%")
println("Quantum speedup: $(round(target_prob * N, digits=2))x")
