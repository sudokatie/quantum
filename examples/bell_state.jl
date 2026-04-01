#!/usr/bin/env julia
# Bell State Example
# Creates all four Bell states and demonstrates entanglement

using Quantum

println("=== Bell States ===\n")

# Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
println("Creating |Φ+⟩ = (|00⟩ + |11⟩)/√2")
c = Circuit(2)
add!(c, H, 1)
add!(c, CNOT, 1, 2)

println("\nCircuit:")
show_circuit(c)

state = execute(c)
println("\nState vector:")
show_state(state)

println("\nProbabilities:")
show_probabilities(state)

println("\nConcurrence (entanglement measure): ", round(concurrence(state), digits=4))
println("Is separable? ", is_separable(state))

# Sample measurements
println("\nSampling 1000 measurements:")
counts = simulate(c, 1000)
for (k, v) in sort(collect(counts))
    basis = lpad(string(k, base=2), 2, '0')
    println("  |$basis⟩: $v times")
end

println("\n=== All Four Bell States ===\n")

# |Φ+⟩ = (|00⟩ + |11⟩)/√2
c1 = Circuit(2)
add!(c1, H, 1)
add!(c1, CNOT, 1, 2)
println("|Φ+⟩:")
show_state(execute(c1))

# |Φ-⟩ = (|00⟩ - |11⟩)/√2
c2 = Circuit(2)
add!(c2, H, 1)
add!(c2, CNOT, 1, 2)
add!(c2, Z, 1)
println("\n|Φ-⟩:")
show_state(execute(c2))

# |Ψ+⟩ = (|01⟩ + |10⟩)/√2
c3 = Circuit(2)
add!(c3, H, 1)
add!(c3, CNOT, 1, 2)
add!(c3, X, 2)
println("\n|Ψ+⟩:")
show_state(execute(c3))

# |Ψ-⟩ = (|01⟩ - |10⟩)/√2
c4 = Circuit(2)
add!(c4, H, 1)
add!(c4, CNOT, 1, 2)
add!(c4, X, 2)
add!(c4, Z, 1)
println("\n|Ψ-⟩:")
show_state(execute(c4))
