# Circuit type for composing quantum operations

"""
    Operation

A single operation in a quantum circuit.
"""
struct Operation
    gate::Matrix{ComplexF64}
    targets::Vector{Int}
    controls::Vector{Int}
    name::String
end

"""
    Circuit

A quantum circuit consisting of a sequence of operations.
"""
mutable struct Circuit
    n_qubits::Int
    operations::Vector{Operation}
end

"""
    Circuit(n::Int) -> Circuit

Create an empty circuit with n qubits.
"""
Circuit(n::Int) = Circuit(n, Operation[])

"""
    add!(circuit::Circuit, gate::Matrix{ComplexF64}, targets::Int...; name="")

Add a gate to the circuit acting on the specified target qubits.
"""
function add!(circuit::Circuit, gate::Matrix{ComplexF64}, targets::Int...; name::String="")
    targets_vec = collect(targets)
    n_gate_qubits = round(Int, log2(size(gate, 1)))
    
    @assert length(targets_vec) == n_gate_qubits "Number of targets must match gate size"
    @assert all(1 .<= targets_vec .<= circuit.n_qubits) "Target qubits out of range"
    
    # Infer name if not provided
    gate_name = name
    if isempty(gate_name)
        gate_name = infer_gate_name(gate)
    end
    
    op = Operation(gate, targets_vec, Int[], gate_name)
    push!(circuit.operations, op)
    circuit
end

"""
    add_controlled!(circuit::Circuit, gate::Matrix{ComplexF64}, 
                    controls::Vector{Int}, targets::Vector{Int}; name="")

Add a controlled gate to the circuit.
"""
function add_controlled!(circuit::Circuit, gate::Matrix{ComplexF64},
                         controls::Vector{Int}, targets::Vector{Int}; name::String="")
    @assert size(gate) == (2, 2) "Controlled gates must use single-qubit base gate"
    @assert all(1 .<= controls .<= circuit.n_qubits) "Control qubits out of range"
    @assert all(1 .<= targets .<= circuit.n_qubits) "Target qubits out of range"
    @assert isempty(intersect(controls, targets)) "Controls and targets must not overlap"
    
    gate_name = name
    if isempty(gate_name)
        base_name = infer_gate_name(gate)
        gate_name = "C" * base_name
    end
    
    op = Operation(gate, targets, controls, gate_name)
    push!(circuit.operations, op)
    circuit
end

"""
    depth(circuit::Circuit) -> Int

Calculate the depth (number of time steps) of the circuit.
Operations on non-overlapping qubits can run in parallel.
"""
function depth(circuit::Circuit)
    if isempty(circuit.operations)
        return 0
    end
    
    # Track the latest time step each qubit is used
    qubit_time = zeros(Int, circuit.n_qubits)
    
    for op in circuit.operations
        all_qubits = vcat(op.targets, op.controls)
        # This operation starts after all its qubits are free
        start_time = maximum(qubit_time[q] for q in all_qubits)
        end_time = start_time + 1
        
        # Update qubit times
        for q in all_qubits
            qubit_time[q] = end_time
        end
    end
    
    maximum(qubit_time)
end

"""
    gate_count(circuit::Circuit) -> Int

Return the total number of gates in the circuit.
"""
gate_count(circuit::Circuit) = length(circuit.operations)

"""
    qubits(circuit::Circuit) -> Int

Return the number of qubits in the circuit.
"""
qubits(circuit::Circuit) = circuit.n_qubits

"""
Infer gate name from matrix.
"""
function infer_gate_name(gate::Matrix{ComplexF64})
    if size(gate) == (2, 2)
        if isapprox(gate, X, atol=1e-10)
            return "X"
        elseif isapprox(gate, Y, atol=1e-10)
            return "Y"
        elseif isapprox(gate, Z, atol=1e-10)
            return "Z"
        elseif isapprox(gate, H, atol=1e-10)
            return "H"
        elseif isapprox(gate, S, atol=1e-10)
            return "S"
        elseif isapprox(gate, T, atol=1e-10)
            return "T"
        elseif isapprox(gate, I_GATE, atol=1e-10)
            return "I"
        else
            return "U"
        end
    elseif size(gate) == (4, 4)
        if isapprox(gate, CNOT, atol=1e-10)
            return "CNOT"
        elseif isapprox(gate, CZ, atol=1e-10)
            return "CZ"
        elseif isapprox(gate, SWAP, atol=1e-10)
            return "SWAP"
        else
            return "U2"
        end
    else
        return "U"
    end
end

# Convenience methods for common gates
add_h!(c::Circuit, t::Int) = add!(c, H, t; name="H")
add_x!(c::Circuit, t::Int) = add!(c, X, t; name="X")
add_y!(c::Circuit, t::Int) = add!(c, Y, t; name="Y")
add_z!(c::Circuit, t::Int) = add!(c, Z, t; name="Z")
add_s!(c::Circuit, t::Int) = add!(c, S, t; name="S")
add_t!(c::Circuit, t::Int) = add!(c, T, t; name="T")
add_cnot!(c::Circuit, ctrl::Int, tgt::Int) = add!(c, CNOT, ctrl, tgt; name="CNOT")
add_cz!(c::Circuit, ctrl::Int, tgt::Int) = add!(c, CZ, ctrl, tgt; name="CZ")
add_swap!(c::Circuit, q1::Int, q2::Int) = add!(c, SWAP, q1, q2; name="SWAP")
