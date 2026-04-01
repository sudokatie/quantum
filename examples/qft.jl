#!/usr/bin/env julia
# Quantum Fourier Transform
# Implements the quantum analog of discrete Fourier transform

using Quantum

println("=== Quantum Fourier Transform ===\n")

# QFT on n qubits
n = 3
println("QFT on $n qubits ($(2^n) basis states)\n")

# Build QFT circuit
function qft_circuit(n::Int)
    c = Circuit(n)
    
    for j in 1:n
        # Hadamard on qubit j
        add!(c, H, j)
        
        # Controlled rotations
        for k in (j+1):n
            # Controlled R_k gate where R_k = diag(1, e^(2πi/2^k))
            angle = 2π / 2^(k - j + 1)
            r_gate = Rz(angle)
            add_controlled!(c, r_gate, [k], [j]; name="R$(k-j+1)")
        end
    end
    
    # Swap qubits to reverse order (standard QFT convention)
    for i in 1:(n÷2)
        add!(c, SWAP, i, n - i + 1)
    end
    
    c
end

c = qft_circuit(n)

println("QFT Circuit:")
show_circuit(c)

println("\nGate count: ", gate_count(c))
println("Circuit depth: ", depth(c))

# Apply QFT to computational basis states
println("\n--- QFT of Basis States ---\n")

for k in 0:min(3, 2^n - 1)
    input = basis_state(n, k)
    output = execute(c, input)
    
    basis = lpad(string(k, base=2), n, '0')
    println("QFT|$basis⟩:")
    
    # Show first few amplitudes
    for i in 0:min(7, 2^n - 1)
        amp = output[i + 1]
        if abs(amp) > 1e-10
            out_basis = lpad(string(i, base=2), n, '0')
            re = round(real(amp), digits=3)
            im = round(imag(amp), digits=3)
            if abs(im) < 1e-10
                println("  |$out_basis⟩: $re")
            else
                sign = im >= 0 ? "+" : "-"
                println("  |$out_basis⟩: $re$sign$(abs(im))i")
            end
        end
    end
    println()
end

# QFT creates uniform superposition from |0...0⟩
println("--- QFT Properties ---\n")

input = zero_state(n)
output = execute(c, input)

println("QFT|000⟩ probabilities:")
show_probabilities(output)

# Verify uniformity
probs = probabilities(output)
expected = 1.0 / 2^n
println("\nExpected uniform probability: $(round(expected, digits=4))")
println("All probabilities equal? ", all(isapprox.(probs, expected, atol=1e-10)))
