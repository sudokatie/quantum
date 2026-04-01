# Quantum gates

"""
Identity gate.
"""
const I_GATE = ComplexF64[1 0; 0 1]

"""
Pauli-X gate (NOT gate, bit flip).
"""
const X = ComplexF64[0 1; 1 0]

"""
Pauli-Y gate.
"""
const Y = ComplexF64[0 -im; im 0]

"""
Pauli-Z gate (phase flip).
"""
const Z = ComplexF64[1 0; 0 -1]

"""
Hadamard gate - creates superposition.
"""
const H = ComplexF64[1 1; 1 -1] / sqrt(2)

"""
Phase gate (S gate, sqrt(Z)).
"""
const S = ComplexF64[1 0; 0 im]

"""
T gate (fourth root of Z).
"""
const T = ComplexF64[1 0; 0 exp(im * pi / 4)]

"""
    Rx(theta) -> Matrix{ComplexF64}

Rotation around X-axis by angle theta.
"""
function Rx(theta::Real)
    c = cos(theta / 2)
    s = sin(theta / 2)
    ComplexF64[c -im*s; -im*s c]
end

"""
    Ry(theta) -> Matrix{ComplexF64}

Rotation around Y-axis by angle theta.
"""
function Ry(theta::Real)
    c = cos(theta / 2)
    s = sin(theta / 2)
    ComplexF64[c -s; s c]
end

"""
    Rz(theta) -> Matrix{ComplexF64}

Rotation around Z-axis by angle theta.
"""
function Rz(theta::Real)
    ComplexF64[exp(-im * theta / 2) 0; 0 exp(im * theta / 2)]
end

"""
    is_unitary(m::Matrix; tol=1e-10) -> Bool

Check if matrix is unitary (M†M = I).
"""
function is_unitary(m::Matrix{ComplexF64}; tol=1e-10)
    n = size(m, 1)
    product = m' * m
    for i in 1:n
        for j in 1:n
            expected = (i == j) ? 1.0 : 0.0
            if abs(product[i, j] - expected) > tol
                return false
            end
        end
    end
    true
end

# Two-qubit gates (Task 4 will expand)

"""
CNOT (Controlled-NOT) gate.
Control on first qubit, target on second.
|00> -> |00>, |01> -> |01>, |10> -> |11>, |11> -> |10>
"""
const CNOT = ComplexF64[
    1 0 0 0
    0 1 0 0
    0 0 0 1
    0 0 1 0
]

"""
CZ (Controlled-Z) gate.
|00> -> |00>, |01> -> |01>, |10> -> |10>, |11> -> -|11>
"""
const CZ = ComplexF64[
    1 0 0 0
    0 1 0 0
    0 0 1 0
    0 0 0 -1
]

"""
SWAP gate - exchanges two qubits.
|00> -> |00>, |01> -> |10>, |10> -> |01>, |11> -> |11>
"""
const SWAP = ComplexF64[
    1 0 0 0
    0 0 1 0
    0 1 0 0
    0 0 0 1
]

"""
iSWAP gate.
"""
const iSWAP = ComplexF64[
    1 0 0 0
    0 0 im 0
    0 im 0 0
    0 0 0 1
]

"""
    tensor_product(a::Matrix, b::Matrix) -> Matrix

Compute Kronecker product of two matrices.
"""
function tensor_product(a::Matrix{ComplexF64}, b::Matrix{ComplexF64})
    kron(a, b)
end

"""
    controlled(gate::Matrix) -> Matrix

Create a controlled version of a single-qubit gate.
Returns a 4x4 matrix with control on first qubit.
"""
function controlled(gate::Matrix{ComplexF64})
    @assert size(gate) == (2, 2) "Gate must be 2x2"
    result = ComplexF64[
        1 0 0 0
        0 1 0 0
        0 0 gate[1,1] gate[1,2]
        0 0 gate[2,1] gate[2,2]
    ]
    result
end
