# optimizer.jl - Variational quantum optimization algorithms
#
# Implements VQE (Variational Quantum Eigensolver) and QAOA
# (Quantum Approximate Optimization Algorithm) for hybrid
# quantum-classical optimization.

using LinearAlgebra

"""
    Hamiltonian

Represents a quantum Hamiltonian as a sum of Pauli terms.
H = Σ_i c_i P_i where P_i are Pauli strings.
"""
struct Hamiltonian
    terms::Vector{Tuple{ComplexF64, Vector{Tuple{Int, Char}}}}  # (coeff, [(qubit, pauli)...])
    n_qubits::Int
end

"""
    hamiltonian(n_qubits::Int, terms::Vector) -> Hamiltonian

Create a Hamiltonian from Pauli terms.
Each term is (coefficient, [(qubit, pauli_char)...])
Pauli chars: 'I', 'X', 'Y', 'Z'
"""
function hamiltonian(n_qubits::Int, terms::Vector)
    Hamiltonian([(ComplexF64(c), collect(ps)) for (c, ps) in terms], n_qubits)
end

"""
    pauli_expectation(state::StateVector, pauli_string::Vector{Tuple{Int, Char}}) -> Float64

Compute expectation value of a Pauli string.
"""
function pauli_expectation(state::StateVector, pauli_string::Vector{Tuple{Int, Char}})
    n = num_qubits(state)
    result = copy(state)
    
    # Apply Pauli operators
    for (qubit, pauli) in pauli_string
        if pauli == 'X'
            apply!(result, X, qubit)
        elseif pauli == 'Y'
            apply!(result, Y, qubit)
        elseif pauli == 'Z'
            apply!(result, Z, qubit)
        end
        # 'I' does nothing
    end
    
    # Compute <ψ|P|ψ>
    real(sum(conj(state.amplitudes) .* result.amplitudes))
end

"""
    expectation_value(state::StateVector, H::Hamiltonian) -> Float64

Compute expectation value <ψ|H|ψ>.
"""
function expectation_value(state::StateVector, H::Hamiltonian)
    total = 0.0
    for (coeff, pauli_string) in H.terms
        exp_val = pauli_expectation(state, pauli_string)
        total += real(coeff) * exp_val
    end
    total
end

"""
    VariationalCircuit

A parameterized quantum circuit for variational algorithms.
"""
struct VariationalCircuit
    n_qubits::Int
    layers::Int
    entanglement::Symbol  # :linear, :full, :circular
end

"""
    variational_circuit(n_qubits::Int; layers=1, entanglement=:linear)

Create a variational ansatz circuit.
"""
function variational_circuit(n_qubits::Int; layers::Int=1, entanglement::Symbol=:linear)
    VariationalCircuit(n_qubits, layers, entanglement)
end

"""
    num_parameters(vc::VariationalCircuit) -> Int

Number of parameters in the variational circuit.
Each layer has 3 rotation parameters per qubit (Rx, Ry, Rz).
"""
function num_parameters(vc::VariationalCircuit)
    vc.layers * vc.n_qubits * 3
end

"""
    build_ansatz(vc::VariationalCircuit, params::Vector{Float64}) -> Circuit

Build a circuit with the given parameters.
"""
function build_ansatz(vc::VariationalCircuit, params::Vector{Float64})
    @assert length(params) == num_parameters(vc) "Expected $(num_parameters(vc)) parameters, got $(length(params))"
    
    circuit = Circuit(vc.n_qubits)
    param_idx = 1
    
    for layer in 1:vc.layers
        # Single-qubit rotations
        for q in 1:vc.n_qubits
            add!(circuit, Rx(params[param_idx]), q)
            add!(circuit, Ry(params[param_idx + 1]), q)
            add!(circuit, Rz(params[param_idx + 2]), q)
            param_idx += 3
        end
        
        # Entangling layer
        if vc.entanglement == :linear
            for q in 1:(vc.n_qubits - 1)
                add!(circuit, CNOT, q, q + 1)
            end
        elseif vc.entanglement == :full
            for q1 in 1:vc.n_qubits
                for q2 in (q1 + 1):vc.n_qubits
                    add!(circuit, CNOT, q1, q2)
                end
            end
        elseif vc.entanglement == :circular
            for q in 1:(vc.n_qubits - 1)
                add!(circuit, CNOT, q, q + 1)
            end
            if vc.n_qubits > 2
                add!(circuit, CNOT, vc.n_qubits, 1)
            end
        end
    end
    
    circuit
end

"""
    OptimizationResult

Result of a variational optimization.
"""
struct OptimizationResult
    optimal_params::Vector{Float64}
    optimal_value::Float64
    iterations::Int
    converged::Bool
    history::Vector{Float64}
end

"""
    nelder_mead(f, x0; max_iter=1000, tol=1e-6)

Simple Nelder-Mead (simplex) optimizer.
"""
function nelder_mead(f, x0::Vector{Float64}; max_iter::Int=1000, tol::Float64=1e-6)
    n = length(x0)
    
    # Initialize simplex
    simplex = [copy(x0)]
    for i in 1:n
        point = copy(x0)
        point[i] += 0.5
        push!(simplex, point)
    end
    
    values = [f(p) for p in simplex]
    history = Float64[minimum(values)]
    
    α, γ, ρ, σ = 1.0, 2.0, 0.5, 0.5  # Standard parameters
    
    for iter in 1:max_iter
        # Sort by function value
        perm = sortperm(values)
        simplex = simplex[perm]
        values = values[perm]
        
        # Check convergence
        if maximum(values) - minimum(values) < tol
            return OptimizationResult(simplex[1], values[1], iter, true, history)
        end
        
        # Centroid (excluding worst)
        centroid = sum(simplex[1:n]) / n
        
        # Reflection
        xr = centroid + α * (centroid - simplex[end])
        fr = f(xr)
        
        if values[1] <= fr < values[n]
            simplex[end] = xr
            values[end] = fr
        elseif fr < values[1]
            # Expansion
            xe = centroid + γ * (xr - centroid)
            fe = f(xe)
            if fe < fr
                simplex[end] = xe
                values[end] = fe
            else
                simplex[end] = xr
                values[end] = fr
            end
        else
            # Contraction
            xc = centroid + ρ * (simplex[end] - centroid)
            fc = f(xc)
            if fc < values[end]
                simplex[end] = xc
                values[end] = fc
            else
                # Shrink
                for i in 2:(n + 1)
                    simplex[i] = simplex[1] + σ * (simplex[i] - simplex[1])
                    values[i] = f(simplex[i])
                end
            end
        end
        
        push!(history, minimum(values))
    end
    
    perm = sortperm(values)
    OptimizationResult(simplex[perm[1]], values[perm[1]], max_iter, false, history)
end

"""
    vqe(H::Hamiltonian, ansatz::VariationalCircuit; 
        initial_params=nothing, max_iter=1000, tol=1e-6) -> OptimizationResult

Variational Quantum Eigensolver.
Finds the ground state energy of Hamiltonian H using variational optimization.
"""
function vqe(H::Hamiltonian, ansatz::VariationalCircuit;
             initial_params::Union{Nothing, Vector{Float64}}=nothing,
             max_iter::Int=1000, tol::Float64=1e-6)
    
    @assert ansatz.n_qubits == H.n_qubits "Ansatz and Hamiltonian must have same number of qubits"
    
    # Initialize parameters
    n_params = num_parameters(ansatz)
    params = isnothing(initial_params) ? 2π * rand(n_params) : copy(initial_params)
    
    # Cost function: expectation value of H
    function cost(p)
        circuit = build_ansatz(ansatz, p)
        state = execute(circuit)
        expectation_value(state, H)
    end
    
    # Optimize
    nelder_mead(cost, params; max_iter=max_iter, tol=tol)
end

"""
    QAOACircuit

QAOA ansatz for optimization problems.
"""
struct QAOACircuit
    n_qubits::Int
    cost_hamiltonian::Hamiltonian
    mixer_hamiltonian::Hamiltonian
    p::Int  # number of layers
end

"""
    maxcut_hamiltonian(edges::Vector{Tuple{Int,Int}}, n_vertices::Int) -> Hamiltonian

Create cost Hamiltonian for MaxCut problem.
H = Σ_{(i,j) ∈ edges} (1 - Z_i Z_j) / 2
"""
function maxcut_hamiltonian(edges::Vector{Tuple{Int,Int}}, n_vertices::Int)
    terms = Tuple{ComplexF64, Vector{Tuple{Int, Char}}}[]
    
    for (i, j) in edges
        # (1 - Z_i Z_j) / 2 = 0.5 * I - 0.5 * Z_i Z_j
        push!(terms, (0.5, Tuple{Int,Char}[]))  # constant term (identity)
        push!(terms, (-0.5, [(i, 'Z'), (j, 'Z')]))
    end
    
    Hamiltonian(terms, n_vertices)
end

"""
    standard_mixer(n_qubits::Int) -> Hamiltonian

Standard QAOA mixer: Σ_i X_i
"""
function standard_mixer(n_qubits::Int)
    terms = [(1.0, [(i, 'X')]) for i in 1:n_qubits]
    Hamiltonian(terms, n_qubits)
end

"""
    qaoa_circuit(cost_H::Hamiltonian, p::Int) -> QAOACircuit

Create QAOA circuit with p layers.
"""
function qaoa_circuit(cost_H::Hamiltonian, p::Int)
    mixer_H = standard_mixer(cost_H.n_qubits)
    QAOACircuit(cost_H.n_qubits, cost_H, mixer_H, p)
end

"""
    build_qaoa(qc::QAOACircuit, gamma::Vector{Float64}, beta::Vector{Float64}) -> Circuit

Build QAOA circuit with given parameters.
gamma: cost Hamiltonian angles
beta: mixer angles
"""
function build_qaoa(qc::QAOACircuit, gamma::Vector{Float64}, beta::Vector{Float64})
    @assert length(gamma) == qc.p "Need $(qc.p) gamma parameters"
    @assert length(beta) == qc.p "Need $(qc.p) beta parameters"
    
    circuit = Circuit(qc.n_qubits)
    
    # Initial superposition
    for q in 1:qc.n_qubits
        add!(circuit, H, q)
    end
    
    # QAOA layers
    for layer in 1:qc.p
        # Cost layer: exp(-i * gamma * H_C)
        # For ZZ terms: exp(-i * gamma * Z_i Z_j) = CNOT; Rz(2*gamma); CNOT
        for (coeff, pauli_string) in qc.cost_hamiltonian.terms
            if length(pauli_string) == 2 && all(p[2] == 'Z' for p in pauli_string)
                q1, q2 = pauli_string[1][1], pauli_string[2][1]
                add!(circuit, CNOT, q1, q2)
                add!(circuit, Rz(2 * gamma[layer] * real(coeff)), q2)
                add!(circuit, CNOT, q1, q2)
            elseif length(pauli_string) == 1 && pauli_string[1][2] == 'Z'
                q = pauli_string[1][1]
                add!(circuit, Rz(2 * gamma[layer] * real(coeff)), q)
            end
            # Identity terms don't affect the state (just global phase)
        end
        
        # Mixer layer: exp(-i * beta * H_M) = Π_i Rx(2*beta)
        for q in 1:qc.n_qubits
            add!(circuit, Rx(2 * beta[layer]), q)
        end
    end
    
    circuit
end

"""
    qaoa(cost_H::Hamiltonian; p=1, max_iter=500, tol=1e-5) -> OptimizationResult

Run QAOA optimization.
"""
function qaoa(cost_H::Hamiltonian; p::Int=1, max_iter::Int=500, tol::Float64=1e-5)
    qc = qaoa_circuit(cost_H, p)
    
    # Cost function
    function cost(params)
        gamma = params[1:p]
        beta = params[p+1:2*p]
        circuit = build_qaoa(qc, gamma, beta)
        state = execute(circuit)
        expectation_value(state, cost_H)
    end
    
    # Initialize randomly
    initial = 2π * rand(2 * p)
    
    nelder_mead(cost, initial; max_iter=max_iter, tol=tol)
end

"""
    qaoa_maxcut(edges::Vector{Tuple{Int,Int}}, n_vertices::Int; p=1, max_iter=500) -> OptimizationResult

Solve MaxCut using QAOA.
"""
function qaoa_maxcut(edges::Vector{Tuple{Int,Int}}, n_vertices::Int; p::Int=1, max_iter::Int=500)
    cost_H = maxcut_hamiltonian(edges, n_vertices)
    result = qaoa(cost_H; p=p, max_iter=max_iter)
    
    # Negate because MaxCut maximizes, but we minimized
    OptimizationResult(
        result.optimal_params,
        -result.optimal_value + 0.5 * length(edges),  # Convert to cut value
        result.iterations,
        result.converged,
        [-v + 0.5 * length(edges) for v in result.history]
    )
end

"""
    sample_qaoa_solution(qc::QAOACircuit, gamma::Vector{Float64}, beta::Vector{Float64}; shots=1000) -> Vector{Int}

Sample bitstrings from the QAOA output state.
Returns the most common bitstring.
"""
function sample_qaoa_solution(qc::QAOACircuit, gamma::Vector{Float64}, beta::Vector{Float64}; shots::Int=1000)
    circuit = build_qaoa(qc, gamma, beta)
    state = execute(circuit)
    
    counts = Dict{Int, Int}()
    for _ in 1:shots
        outcome = sample(state)
        counts[outcome] = get(counts, outcome, 0) + 1
    end
    
    # Return most common
    best = argmax(counts)
    best
end
