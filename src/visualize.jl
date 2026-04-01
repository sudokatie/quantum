# Visualization functions

"""
    show_circuit(io::IO, circuit::Circuit)

Display an ASCII circuit diagram.
"""
function show_circuit(io::IO, circuit::Circuit)
    if isempty(circuit.operations)
        println(io, "Empty circuit with $(circuit.n_qubits) qubits")
        return
    end
    
    n = circuit.n_qubits
    
    # Build the circuit diagram
    # Each qubit gets a line, operations are columns
    lines = [String[] for _ in 1:n]
    
    # Track which column each qubit is at
    qubit_col = ones(Int, n)
    
    for op in circuit.operations
        all_qubits = sort(vcat(op.targets, op.controls))
        max_col = maximum(qubit_col[q] for q in all_qubits)
        
        # Pad other involved qubits to same column
        for q in all_qubits
            while qubit_col[q] < max_col
                push!(lines[q], "───")
                qubit_col[q] += 1
            end
        end
        
        # Add the operation
        if length(op.targets) == 1 && isempty(op.controls)
            # Single qubit gate
            t = op.targets[1]
            name = length(op.name) == 1 ? "[$(op.name)]" : 
                   length(op.name) == 2 ? "[$(op.name)]" : "[$(op.name[1:min(2,length(op.name))])]"
            push!(lines[t], name)
            qubit_col[t] += 1
        elseif !isempty(op.controls)
            # Controlled gate
            min_q = minimum(all_qubits)
            max_q = maximum(all_qubits)
            
            for q in min_q:max_q
                if q in op.controls
                    push!(lines[q], "─●─")
                elseif q in op.targets
                    name = op.gate == X || isapprox(op.gate, X, atol=1e-10) ? "─⊕─" : "─[$(op.name[1])]─"
                    push!(lines[q], name)
                else
                    push!(lines[q], "─│─")
                end
                qubit_col[q] += 1
            end
        elseif length(op.targets) == 2
            # Two-qubit gate (CNOT, SWAP, etc.)
            t1, t2 = op.targets
            min_q = min(t1, t2)
            max_q = max(t1, t2)
            
            for q in min_q:max_q
                if op.name == "CNOT"
                    if q == t1
                        push!(lines[q], "─●─")
                    elseif q == t2
                        push!(lines[q], "─⊕─")
                    else
                        push!(lines[q], "─│─")
                    end
                elseif op.name == "SWAP"
                    if q == t1 || q == t2
                        push!(lines[q], "─╳─")
                    else
                        push!(lines[q], "─│─")
                    end
                else
                    if q == t1 || q == t2
                        push!(lines[q], "─■─")
                    else
                        push!(lines[q], "─│─")
                    end
                end
                qubit_col[q] += 1
            end
        end
    end
    
    # Pad all lines to same length
    max_len = maximum(length(l) for l in lines)
    for i in 1:n
        while length(lines[i]) < max_len
            push!(lines[i], "───")
        end
    end
    
    # Print
    for (i, line) in enumerate(lines)
        print(io, "q$(i): ")
        println(io, join(line, ""))
    end
end

show_circuit(circuit::Circuit) = show_circuit(stdout, circuit)

"""
    show_state(io::IO, state::StateVector; threshold=1e-10)

Display state vector in ket notation, showing only non-zero amplitudes.
"""
function show_state(io::IO, state::StateVector; threshold::Float64=1e-10)
    n = state.n_qubits
    
    for i in 0:(2^n - 1)
        amp = state.amplitudes[i + 1]
        if abs(amp) > threshold
            # Format amplitude
            re = real(amp)
            im = imag(amp)
            
            if abs(im) < threshold
                amp_str = @sprintf("%.4f", re)
            elseif abs(re) < threshold
                amp_str = @sprintf("%.4fi", im)
            else
                sign = im >= 0 ? "+" : "-"
                amp_str = @sprintf("%.4f%s%.4fi", re, sign, abs(im))
            end
            
            # Format basis state
            basis = lpad(string(i, base=2), n, '0')
            println(io, "|$(basis)⟩: $(amp_str)")
        end
    end
end

show_state(state::StateVector; threshold::Float64=1e-10) = 
    show_state(stdout, state; threshold=threshold)

"""
    show_probabilities(io::IO, state::StateVector; threshold=1e-10, bar_width=40)

Display probability distribution as histogram.
"""
function show_probabilities(io::IO, state::StateVector; 
                            threshold::Float64=1e-10, bar_width::Int=40)
    n = state.n_qubits
    probs = probabilities(state)
    
    for i in 0:(2^n - 1)
        p = probs[i + 1]
        if p > threshold
            basis = lpad(string(i, base=2), n, '0')
            bar_len = round(Int, p * bar_width)
            bar = repeat("█", bar_len)
            println(io, "|$(basis)⟩: $(bar) $(round(p * 100, digits=1))%")
        end
    end
end

show_probabilities(state::StateVector; threshold::Float64=1e-10, bar_width::Int=40) = 
    show_probabilities(stdout, state; threshold=threshold, bar_width=bar_width)

"""
    show_bloch(io::IO, state::StateVector)

Display Bloch sphere representation for a single qubit state.
"""
function show_bloch(io::IO, state::StateVector)
    @assert state.n_qubits == 1 "Bloch sphere only for single qubit"
    
    # Calculate Bloch vector components
    # For |psi> = a|0> + b|1>
    # x = 2*Re(a*conj(b))
    # y = 2*Im(a*conj(b))
    # z = |a|^2 - |b|^2
    
    a = state.amplitudes[1]
    b = state.amplitudes[2]
    
    x = 2 * real(a * conj(b))
    y = 2 * imag(a * conj(b))
    z = abs2(a) - abs2(b)
    
    println(io, "Bloch vector:")
    println(io, "  x = $(round(x, digits=4))")
    println(io, "  y = $(round(y, digits=4))")
    println(io, "  z = $(round(z, digits=4))")
    
    # Simple ASCII art
    println(io, "")
    println(io, "      ┌───────┐")
    z_pos = z > 0.5 ? "  ●   " : z > 0 ? "   ●  " : z > -0.5 ? "    ● " : "     ●"
    println(io, "  |0⟩ │$(z_pos)│")
    println(io, "      │       │")
    println(io, "      │       │")
    println(io, "  |1⟩ └───────┘")
end

show_bloch(state::StateVector) = show_bloch(stdout, state)
