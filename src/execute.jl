# Circuit execution

"""
    execute(circuit::Circuit) -> StateVector

Execute the circuit starting from |0>^n.
"""
function execute(circuit::Circuit)
    execute(circuit, zero_state(circuit.n_qubits))
end

"""
    execute(circuit::Circuit, initial::StateVector) -> StateVector

Execute the circuit starting from the given initial state.
"""
function execute(circuit::Circuit, initial::StateVector)
    @assert initial.n_qubits == circuit.n_qubits "State size must match circuit"
    
    state = initial
    
    for op in circuit.operations
        if isempty(op.controls)
            # Regular gate application
            if length(op.targets) == 1
                state = apply(state, op.gate, op.targets[1])
            else
                state = apply(state, op.gate, op.targets)
            end
        else
            # Controlled gate - need to build the full controlled matrix
            state = apply_controlled_op(state, op)
        end
    end
    
    state
end

"""
Apply a controlled operation to a state.
"""
function apply_controlled_op(state::StateVector, op::Operation)
    n = state.n_qubits
    dim = 2^n
    new_amps = copy(state.amplitudes)
    
    # For each basis state, check if all controls are 1
    # If so, apply the gate to target
    for i in 0:(dim - 1)
        # Check all controls are 1
        all_controls_set = true
        for c in op.controls
            bit_pos = n - c
            if (i >> bit_pos) & 1 == 0
                all_controls_set = false
                break
            end
        end
        
        if !all_controls_set
            continue
        end
        
        # Apply single-qubit gate to target
        # Only process when target bit is 0 to avoid double processing
        for t in op.targets
            bit_pos = n - t
            if (i >> bit_pos) & 1 == 0
                j = i | (1 << bit_pos)
                
                a0 = state.amplitudes[i + 1]
                a1 = state.amplitudes[j + 1]
                
                new_amps[i + 1] = op.gate[1,1] * a0 + op.gate[1,2] * a1
                new_amps[j + 1] = op.gate[2,1] * a0 + op.gate[2,2] * a1
            end
        end
    end
    
    StateVector(new_amps)
end

"""
    execute_and_measure(circuit::Circuit; shots::Int=1000) -> Dict{Int, Int}

Execute the circuit and measure, returning measurement counts.
"""
function execute_and_measure(circuit::Circuit; shots::Int=1000)
    state = execute(circuit)
    sample(state, shots)
end

"""
    execute_and_measure(circuit::Circuit, initial::StateVector; shots::Int=1000) -> Dict{Int, Int}

Execute from initial state and measure.
"""
function execute_and_measure(circuit::Circuit, initial::StateVector; shots::Int=1000)
    state = execute(circuit, initial)
    sample(state, shots)
end

"""
    simulate(circuit::Circuit, shots::Int) -> Dict{Int, Int}

Alias for execute_and_measure.
"""
simulate(circuit::Circuit, shots::Int) = execute_and_measure(circuit; shots=shots)

"""
    simulate(circuit::Circuit, initial::StateVector, shots::Int) -> Dict{Int, Int}

Simulate from initial state.
"""
simulate(circuit::Circuit, initial::StateVector, shots::Int) = 
    execute_and_measure(circuit, initial; shots=shots)
