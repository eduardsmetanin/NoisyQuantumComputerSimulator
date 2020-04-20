using Test, LinearAlgebra, NoisyQuantumComputerSimulator, NoisyQuantumComputerSimulator.FullRegisterGate

⊗ = kron

I = Matrix{Complex{Float64}}([1 0; 0 1])
H = Matrix{Complex{Float64}}(1/√2 * [1 1; 1 -1])
X = Matrix{Complex{Float64}}([0 1; 1 0])
S = Matrix{Complex{Float64}}([1 0; 0 complex(0, 1)])

# Testing Program.

p = Program(1)
@test p.density_matrix == [1 0; 0 0]
p = Program(2)
@test p.density_matrix == [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
p = Program(3)
@test p.density_matrix == [1 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0]

# Testing CNOT implemented via x.

BELL = 1/√2 * [1 0 1 0; 0 1 0 1; 0 1 0 -1; 1 0 -1 0] # (CX)(H⊗I)
BELL₀ = BELL * [1; 0; 0; 0]
BELL₁ = BELL * [0; 1; 0; 0]
BELL₂ = BELL * [0; 0; 1; 0]
BELL₃ = BELL * [0; 0; 0; 1]

p = Program(2)
h(p,1)
x(p,0,[1])
exec(p)
@test p.density_matrix == BELL₀ * BELL₀'

p = Program(2)
x(p,0)
h(p,1)
x(p,0,[1])
exec(p)
@test p.density_matrix == BELL₁ * BELL₁'

p = Program(2)
x(p,1)
h(p,1)
x(p,0,[1])
exec(p)
@test p.density_matrix == BELL₂ * BELL₂'

p = Program(2)
x(p,0)
x(p,1)
h(p,1)
x(p,0,[1])
exec(p)
@test p.density_matrix == BELL₃ * BELL₃'

p = Program(4)
h(p,2)
x(p,1,[2])
exec(p)
@test p.density_matrix ≈ ([1; 0] ⊗ BELL₀ ⊗ [1; 0]) * ([1 0] ⊗ BELL₀' ⊗ [1 0])

# Testing CCNOT implemented via x.

CCBELL = [.5  .5   0   0  .5  .5   0   0;
          .5 -.5   0   0  .5 -.5   0   0;
           0   0  .5  .5   0   0  .5  .5;
           0   0  .5 -.5   0   0  .5 -.5;
          .5  .5   0   0 -.5 -.5   0   0;
           0   0  .5 -.5   0   0 -.5  .5;
           0   0  .5  .5   0   0 -.5 -.5;
          .5 -.5   0   0 -.5  .5   0   0]
CCBELL₀ = CCBELL * [1; 0; 0; 0; 0; 0; 0; 0]
CCBELL₁ = CCBELL * [0; 1; 0; 0; 0; 0; 0; 0]
CCBELL₂ = CCBELL * [0; 0; 1; 0; 0; 0; 0; 0]
CCBELL₃ = CCBELL * [0; 0; 0; 1; 0; 0; 0; 0]
CCBELL₄ = CCBELL * [0; 0; 0; 0; 1; 0; 0; 0]
CCBELL₅ = CCBELL * [0; 0; 0; 0; 0; 1; 0; 0]
CCBELL₆ = CCBELL * [0; 0; 0; 0; 0; 0; 1; 0]
CCBELL₇ = CCBELL * [0; 0; 0; 0; 0; 0; 0; 1]

p = Program(3)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₀ * CCBELL₀'

p = Program(3)
x(p,0)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₁ * CCBELL₁'

p = Program(3)
x(p,1)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₂ * CCBELL₂'

p = Program(3)
x(p,0)
x(p,1)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₃ * CCBELL₃'

p = Program(3)
x(p,2)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₄ * CCBELL₄'

p = Program(3)
x(p,0)
x(p,2)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₅ * CCBELL₅'

p = Program(3)
x(p,1)
x(p,2)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₆ * CCBELL₆'

p = Program(3)
x(p,0)
x(p,1)
x(p,2)
h(p,0)
h(p,2)
x(p,1,[0,2])
exec(p)
@test p.density_matrix ≈ CCBELL₇ * CCBELL₇'

# Testing h.

p = Program(1)
h(p,0)
exec(p)
@test p.density_matrix ≈ [.5 .5; .5 .5]

# Testing x.

p = Program(1)
x(p,0)
exec(p)
@test p.density_matrix == [0 0; 0 1]

p = Program(1)
gate(p,X,0)
x(p,0)
exec(p)
@test p.density_matrix == [1 0; 0 0]

p = Program(1)
h(p,0)
x(p,0)
exec(p)
@test p.density_matrix ≈ [.5 .5; .5 .5]

p = Program(1)
gate(p,X,0)
h(p,0)
x(p,0)
exec(p)
@test p.density_matrix ≈ [.5 -.5; -.5 .5]

# Testing y.

p = Program(1)
y(p,0)
exec(p)
@test p.density_matrix == [0 0; 0 1]

p = Program(1)
x(p,0)
y(p,0)
exec(p)
@test p.density_matrix == [1 0; 0 0]

# Testing z.

p = Program(1)
z(p,0)
exec(p)
@test p.density_matrix == [1 0; 0 0]

p = Program(1)
x(p,0)
z(p,0)
exec(p)
@test p.density_matrix == [0 0; 0 1]

p = Program(1)
h(p,0)
z(p,0)
exec(p)
@test p.density_matrix ≈ [.5 -.5; -.5 .5]

p = Program(1)
x(p,0)
h(p,0)
z(p,0)
exec(p)
@test p.density_matrix ≈ [.5 .5; .5 .5]

# Testing swap.

p = Program(2)
x(p,0)
swap(p,0,1)
exec(p)
@test p.density_matrix == [0; 0; 1; 0] * [0 0 1 0]

p = Program(4)
x(p,1)
x(p,0)
x(p,3)
swap(p,1,2,[0,3])
exec(p)
@test p.density_matrix == [0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1] * ([0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1])'

p = Program(4)
x(p,1)
x(p,3)
swap(p,1,2,[0,3]) # Control-SWAP does nothing because one of the control qubits is 0.
exec(p)
@test p.density_matrix == [0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0] * ([0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0])'

p = Program(2)
x(p,0)
h(p,0)
swap(p,0,1)
exec(p)
@test p.density_matrix == [1/√2; 0; -1/√2; 0] * [1/√2 0 -1/√2 0]

p = Program(2)
x(p,0)
h(p,0)
swap(p,0,1)
exec(p)
@test p.density_matrix == [1/√2; 0; -1/√2; 0] * [1/√2 0 -1/√2 0]

p = Program(3)
x(p,0)
swap(p,0,1,[2])
exec(p)
@test p.density_matrix == [0; 1; 0; 0; 0; 0; 0; 0] * [0 1 0 0 0 0 0 0]

# Testing s.

p = Program(1)
h(p,0)
s(p,0)
exec(p)
@test p.density_matrix ≈ [.5 -.5im; .5im .5]

p = Program(1)
x(p,0)
h(p,0)
s(p,0)
exec(p)
@test p.density_matrix ≈ [.5 .5im; -.5im .5]

# Testing amplitude damping. Amplitude damping causes |1> to decay to |0> with some probability.

# When we start in |0> state, damping kraus with any |1> to |0> decay probability never gets applied,
# residual kraus always gets applied, so the state doesn't change.
p = Program(1)
amplitude_damping(p, 1, 0)
exec(p)
@test p.density_matrix == [1 0; 0 0]

# Starting in |1> state, setting |1> to |0> decay probability to 0.
# Residual kraus should be applied and the state should remain |1>.
p = Program(1)
x(p,0)
amplitude_damping(p, 0, 0)
exec(p)
@test p.density_matrix == [0 0; 0 1]

# Starting in |1> state, setting |1> to |0> decay probability to 1.
# Damping kraus should be applied and the state should decay to |0>.
p = Program(1)
x(p,0)
amplitude_damping(p, 1, 0)
exec(p)
@test p.density_matrix == [1 0; 0 0]

# Starting in |+> state, setting |1> to |0> decay probability to 0.
# Residual kraus should be applied and the state should remain |1>.
p = Program(1)
h(p,0)
amplitude_damping(p, 0, 0)
exec(p)
@test p.density_matrix ≈ [.5 .5; .5 .5]

# Starting in |+> state, setting |1> to |0> decay probability to 1.
# Damping kraus should be applied and the state should decay to |0>.
p = Program(1)
h(p,0)
amplitude_damping(p, 1, 0)
exec(p)
@test p.density_matrix ≈ [1 0; 0 0]

# Testing private functions.

CX = [1 0 0 0;
      0 1 0 0;
      0 0 0 1;
      0 0 1 0]

CCX = [1 0 0 0 0 0 0 0;
       0 1 0 0 0 0 0 0;
       0 0 1 0 0 0 0 0;
       0 0 0 1 0 0 0 0;
       0 0 0 0 1 0 0 0;
       0 0 0 0 0 1 0 0;
       0 0 0 0 0 0 0 1;
       0 0 0 0 0 0 1 0]

CCCX = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]

XC = [1 0 0 0;
      0 0 0 1;
      0 0 1 0;
      0 1 0 0]

XCC = [1 0 0 0 0 0 0 0;
       0 1 0 0 0 0 0 0;
       0 0 1 0 0 0 0 0;
       0 0 0 0 0 0 0 1;
       0 0 0 0 1 0 0 0;
       0 0 0 0 0 1 0 0;
       0 0 0 0 0 0 1 0;
       0 0 0 1 0 0 0 0]

CXC = [1 0 0 0 0 0 0 0;
       0 1 0 0 0 0 0 0;
       0 0 1 0 0 0 0 0;
       0 0 0 1 0 0 0 0;
       0 0 0 0 1 0 0 0;
       0 0 0 0 0 0 0 1;
       0 0 0 0 0 0 1 0;
       0 0 0 0 0 1 0 0]

CIXC = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]

CXIC = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]

τ = 1/√2
CHIC = [1 0 0 0 0 0 0 0 0 0 0 0 0  0 0  0;
        0 1 0 0 0 0 0 0 0 0 0 0 0  0 0  0;
        0 0 1 0 0 0 0 0 0 0 0 0 0  0 0  0;
        0 0 0 1 0 0 0 0 0 0 0 0 0  0 0  0;
        0 0 0 0 1 0 0 0 0 0 0 0 0  0 0  0;
        0 0 0 0 0 1 0 0 0 0 0 0 0  0 0  0;
        0 0 0 0 0 0 1 0 0 0 0 0 0  0 0  0;
        0 0 0 0 0 0 0 1 0 0 0 0 0  0 0  0;
        0 0 0 0 0 0 0 0 1 0 0 0 0  0 0  0;
        0 0 0 0 0 0 0 0 0 τ 0 0 0  τ 0  0;
        0 0 0 0 0 0 0 0 0 0 1 0 0  0 0  0;
        0 0 0 0 0 0 0 0 0 0 0 τ 0  0 0  τ;
        0 0 0 0 0 0 0 0 0 0 0 0 1  0 0  0;
        0 0 0 0 0 0 0 0 0 τ 0 0 0 -τ 0  0;
        0 0 0 0 0 0 0 0 0 0 0 0 0  0 1  0;
        0 0 0 0 0 0 0 0 0 0 0 τ 0  0 0 -τ]

HH = Matrix{Complex{Float64}}([.5 .5 .5 .5; .5 -.5 .5 -.5; .5 .5 -.5 -.5; .5 -.5 -.5 .5])

CHH = [1 0 0 0  0   0   0   0;
       0 1 0 0  0   0   0   0;
       0 0 1 0  0   0   0   0;
       0 0 0 1  0   0   0   0;
       0 0 0 0 .5  .5  .5  .5;
       0 0 0 0 .5 -.5  .5 -.5;
       0 0 0 0 .5  .5 -.5 -.5;
       0 0 0 0 .5 -.5 -.5  .5]

@test FullRegisterGate.build(1, H, 0) == H
@test FullRegisterGate.build(2, HH, 0) == HH
@test FullRegisterGate.build(3, HH, 0, [2]) == CHH
@test FullRegisterGate.build(6, HH, 2, [4]) == I ⊗ CHH ⊗ I ⊗ I
@test FullRegisterGate.build(6, H, 3, [4,1]) == I ⊗ CHIC ⊗ I
@test FullRegisterGate.build(2, X, 0, [1]) == CX
@test FullRegisterGate.build(3, X, 0, [1,2]) == CCX
@test FullRegisterGate.build(4, X, 0, [2,1,3]) == CCCX
@test FullRegisterGate.build(3, X, 0, [1]) == I ⊗ CX
@test FullRegisterGate.build(3, X, 1, [2]) == CX ⊗ I
@test FullRegisterGate.build(4, X, 0, [1]) == I ⊗ I ⊗ CX
@test FullRegisterGate.build(5, X, 0, [1,2]) == I ⊗ I ⊗ CCX
@test FullRegisterGate.build(2, X, 1, [0]) == XC
@test FullRegisterGate.build(3, X, 1, [0]) == I ⊗ XC
@test FullRegisterGate.build(3, X, 2, [0,1]) == XCC
@test FullRegisterGate.build(4, X, 2, [1,0]) == I ⊗ XCC
@test FullRegisterGate.build(3, X, 1, [0,2]) == CXC
@test FullRegisterGate.build(4, X, 2, [1,3]) == CXC ⊗ I
@test FullRegisterGate.build(5, X, 2, [3,1]) == I ⊗ CXC ⊗ I
@test FullRegisterGate.build(4, X, 1, [0,3]) == CIXC
@test FullRegisterGate.build(3, X, 2, [1]) == XC ⊗ I
@test FullRegisterGate.build(4, X, 2, [3,0]) == CXIC
@test FullRegisterGate.build(7, X, 5, [3,6]) == CXIC ⊗ I ⊗ I ⊗ I
@test FullRegisterGate.build(6, X, 3, [1,4]) == I ⊗ CXIC ⊗ I
@test FullRegisterGate.build(1, S, 0) == S
#Testing asymmetric gates to ensure the gate doesn't get transposed by mistake.
@test FullRegisterGate.build(1, NoisyQuantumComputerSimulator.damping_residual_kraus(), 0) == NoisyQuantumComputerSimulator.damping_residual_kraus()
@test FullRegisterGate.build(3, NoisyQuantumComputerSimulator.damping_residual_kraus(), 1) == I ⊗ NoisyQuantumComputerSimulator.damping_residual_kraus() ⊗ I
@test FullRegisterGate.build(4, NoisyQuantumComputerSimulator.damping_residual_kraus(), 2) == I ⊗ NoisyQuantumComputerSimulator.damping_residual_kraus() ⊗ I ⊗ I
@test FullRegisterGate.build(1, NoisyQuantumComputerSimulator.damping_kraus(), 0) == NoisyQuantumComputerSimulator.damping_kraus()
@test FullRegisterGate.build(3, NoisyQuantumComputerSimulator.damping_kraus(), 1) == I ⊗ NoisyQuantumComputerSimulator.damping_kraus() ⊗ I
@test FullRegisterGate.build(4, NoisyQuantumComputerSimulator.damping_kraus(), 1) == I ⊗ I ⊗ NoisyQuantumComputerSimulator.damping_kraus() ⊗ I

@test FullRegisterGate.gate_full_register_bitmask(1, 0) == 0b1
@test FullRegisterGate.gate_full_register_bitmask(1, 1) == 0b10
@test FullRegisterGate.gate_full_register_bitmask(1, 2) == 0b100
@test FullRegisterGate.gate_full_register_bitmask(2, 0) == 0b11
@test FullRegisterGate.gate_full_register_bitmask(2, 1) == 0b110
@test FullRegisterGate.gate_full_register_bitmask(2, 2) == 0b1100
@test FullRegisterGate.gate_full_register_bitmask(3, 0) == 0b111
@test FullRegisterGate.gate_full_register_bitmask(3, 1) == 0b1110
@test FullRegisterGate.gate_full_register_bitmask(3, 2) == 0b11100

@test FullRegisterGate.control_full_register_bitmask(Array{Int64,1}([])) == 0
@test FullRegisterGate.control_full_register_bitmask([0]) == 1
@test FullRegisterGate.control_full_register_bitmask([0, 1]) == 0b11
@test FullRegisterGate.control_full_register_bitmask([1, 5]) == 0b100010
@test FullRegisterGate.control_full_register_bitmask([1, 7, 5]) == 0b10100010

@test FullRegisterGate.all_control_bits_set(0b11, 0b10) == true
@test FullRegisterGate.all_control_bits_set(0b1100, 0b1010) == false
@test FullRegisterGate.all_control_bits_set(0b0000, 0b1010) == false
@test FullRegisterGate.all_control_bits_set(0b1111, 0b0101) == true
@test FullRegisterGate.all_control_bits_set(0b0, 0b0) == true
@test FullRegisterGate.all_control_bits_set(0, 0) == true
@test FullRegisterGate.all_control_bits_set(0xffffffffffffffff, 0xffffffffffffffff) == true
@test FullRegisterGate.all_control_bits_set(0x8fffffffffffffff, 0xffffffffffffffff) == false

@test NoisyQuantumComputerSimulator.damping_kraus_operators(.1) == [[1 0; 0 √.9], [0 √.1; 0 0]]

@test NoisyQuantumComputerSimulator.value_index_by_probabilities(.0, [.3, .4, .2, .1]) == 1
@test NoisyQuantumComputerSimulator.value_index_by_probabilities(.35, [.3, .4, .2, .1]) == 2
@test NoisyQuantumComputerSimulator.value_index_by_probabilities(.8, [.3, .4, .2, .1]) == 3
@test NoisyQuantumComputerSimulator.value_index_by_probabilities(.95, [.3, .4, .2, .1]) == 4
@test NoisyQuantumComputerSimulator.value_index_by_probabilities(.9999, [.3, .4, .2, .1]) == 4
@test NoisyQuantumComputerSimulator.value_index_by_probabilities(1.0, [.3, .4, .2, .1]) == 4

# @test product!([1; 0], I) == [1; 0]
# @test product!([0; 1], I) == [0; 1]
# @test product!([1; 0], H) == [1/√2; 1/√2]
# @test product!([0; 1], H) == [1/√2; -1/√2]
# @test product!([1 0; 0 0], I) == [1 0; 0 0]
# @test product!([0 0; 0 1], I) == [0 0; 0 1]
# @test product!([1 0; 0 0], H) ≈ [0.5 0.5; 0.5 0.5]
# @test product!([0 0; 0 1], H) ≈ [0.5 -0.5; -0.5 0.5]
# @test product!([0.5 0.5; 0.5 0.5], H) ≈ [1 0; 0 0]
# @test product!([0.5 -0.5; -0.5 0.5], H) ≈ [0 0; 0 1]
