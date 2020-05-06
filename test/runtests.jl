using Test
using LinearAlgebra
using NoisyQuantumComputerSimulator
using NoisyQuantumComputerSimulator.FullRegisterGate
using NoisyQuantumComputerSimulator.Gates

⊗ = kron

matrix_I = Matrix{Complex{Float64}}([1 0; 0 1])
matrix_H = Matrix{Complex{Float64}}(1/√2 * [1 1; 1 -1])
matrix_X = Matrix{Complex{Float64}}([0 1; 1 0])
matrix_S = Matrix{Complex{Float64}}([1 0; 0 complex(0, 1)])

# Testing Curcuit.

c = Curcuit(1)
@test c.density_matrix == [1 0; 0 0]
c = Curcuit(2)
@test c.density_matrix == [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
c = Curcuit(3)
@test c.density_matrix == [1 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0;
                           0 0 0 0 0 0 0 0]

# Testing CNOT.

BELL = 1/√2 * [1 0 1 0; 0 1 0 1; 0 1 0 -1; 1 0 -1 0] # (CX)(H⊗I)
BELL₀ = BELL * [1; 0; 0; 0]
BELL₁ = BELL * [0; 1; 0; 0]
BELL₂ = BELL * [0; 0; 1; 0]
BELL₃ = BELL * [0; 0; 0; 1]

c = Curcuit(2)
c += H(1)
c += CNOT(1,0)
exec(c)
@test c.density_matrix == BELL₀ * BELL₀'

c = Curcuit(2)
c += X(0)
c += H(1)
c += CNOT(1,0)
exec(c)
@test c.density_matrix == BELL₁ * BELL₁'

c = Curcuit(2)
c += X(1)
c += H(1)
c += CNOT(1,0)
exec(c)
@test c.density_matrix == BELL₂ * BELL₂'

c = Curcuit(2)
c += X(0)
c += X(1)
c += H(1)
c += CNOT(1,0)
exec(c)
@test c.density_matrix == BELL₃ * BELL₃'

c = Curcuit(4)
c += H(2)
c += CNOT(2,1)
exec(c)
@test c.density_matrix ≈ ([1; 0] ⊗ BELL₀ ⊗ [1; 0]) * ([1 0] ⊗ BELL₀' ⊗ [1 0])

# Testing CCNOT.

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

c = Curcuit(3)
c += H(0)
c += H(2)
c += CCNOT(0,2,1)
exec(c)
@test c.density_matrix ≈ CCBELL₀ * CCBELL₀'

c = Curcuit(3)
c += X(0)
c += H(0)
c += H(2)
c += CCNOT(0,2,1)
exec(c)
@test c.density_matrix ≈ CCBELL₁ * CCBELL₁'

c = Curcuit(3)
c += X(1)
c += H(0)
c += H(2)
c += CCNOT(0,2,1)
exec(c)
@test c.density_matrix ≈ CCBELL₂ * CCBELL₂'

c = Curcuit(3)
c += X(0)
c += X(1)
c += H(0)
c += H(2)
c += CCNOT(0,2,1)
exec(c)
@test c.density_matrix ≈ CCBELL₃ * CCBELL₃'

c = Curcuit(3)
c += X(2)
c += H(0)
c += H(2)
c += CCNOT(0,2,1)
exec(c)
@test c.density_matrix ≈ CCBELL₄ * CCBELL₄'

c = Curcuit(3)
c += X(0)
c += X(2)
c += H(0)
c += H(2)
c += CCNOT(0,2,1)
exec(c)
@test c.density_matrix ≈ CCBELL₅ * CCBELL₅'

c = Curcuit(3)
c += X(1)
c += X(2)
c += H(0)
c += H(2)
c += CCNOT(0,2,1)
exec(c)
@test c.density_matrix ≈ CCBELL₆ * CCBELL₆'

c = Curcuit(3)
c += X(0)
c += X(1)
c += X(2)
c += H(0)
c += H(2)
# c += CCNOT(0,2,1)
c += custom_gate(matrix_X,1,[0,2]) #c += custom_gate(X,0)
exec(c)
@test c.density_matrix ≈ CCBELL₇ * CCBELL₇'

# Testing X.

c = Curcuit(1)
c += X(0)
exec(c)
@test c.density_matrix == [0 0; 0 1]

c = Curcuit(1)
c += custom_gate(matrix_X,0)
c += X(0)
exec(c)
@test c.density_matrix == [1 0; 0 0]

c = Curcuit(1)
c += H(0)
c += X(0)
exec(c)
@test c.density_matrix ≈ [.5 .5; .5 .5]

c = Curcuit(1)
c += custom_gate(matrix_X,0)
c += H(0)
c += X(0)
exec(c)
@test c.density_matrix ≈ [.5 -.5; -.5 .5]

# Testing Y.

c = Curcuit(1)
c += Y(0)
exec(c)
@test c.density_matrix == [0 0; 0 1]

c = Curcuit(1)
c += X(0)
c += Y(0)
exec(c)
@test c.density_matrix == [1 0; 0 0]

# Testing Z.

c = Curcuit(1)
c += Z(0)
exec(c)
@test c.density_matrix == [1 0; 0 0]

c = Curcuit(1)
c += X(0)
c += Z(0)
exec(c)
@test c.density_matrix == [0 0; 0 1]

c = Curcuit(1)
c += H(0)
c += Z(0)
exec(c)
@test c.density_matrix ≈ [.5 -.5; -.5 .5]

c = Curcuit(1)
c += X(0)
c += H(0)
c += Z(0)
exec(c)
@test c.density_matrix ≈ [.5 .5; .5 .5]

# Testing H.

c = Curcuit(1)
c += H(0)
exec(c)
@test c.density_matrix ≈ [.5 .5; .5 .5]

# Testing PHASE.

c = Curcuit(1)
c += H(0)
c += PHASE(pi/2, 0)
exec(c)
@test c.density_matrix ≈ [.5 complex(0,-.5); complex(0,.5) .5]

# Testing S.

c = Curcuit(1)
c += H(0)
c += S(0)
exec(c)
@test c.density_matrix ≈ [.5 -.5im; .5im .5]

c = Curcuit(1)
c += X(0)
c += H(0)
c += S(0)
exec(c)
@test c.density_matrix ≈ [.5 .5im; -.5im .5]

# Testing T.

c = Curcuit(1)
c += H(0)
c += T(0)
exec(c)
@test c.density_matrix ≈ [.5 (1-1im)/(2*sqrt(2)); (1+1im)/(2*sqrt(2)) .5]

# Testing CZ.

c = Curcuit(2)
c += H(0)
c += H(1)
c += CZ(1,0)
exec(c)
@test c.density_matrix ≈ [.5; .5; .5; -.5] * [.5 .5 .5 -.5]

# Testing RX.
c = Curcuit(1)
c += RX(pi,0)
exec(c)
@test c.density_matrix ≈ [0 0; 0 1]

c = Curcuit(1)
c += RX(pi/2,0)
exec(c)
@test c.density_matrix ≈ [.5 .5im; -.5im .5]

# Testing RY.
c = Curcuit(1)
c += RY(pi,0)
exec(c)
@test c.density_matrix ≈ [0 0; 0 1]

c = Curcuit(1)
c += RY(pi/2,0)
exec(c)
@test c.density_matrix ≈ [.5 .5; .5 .5]

# Testing RZ.

c = Curcuit(1)
c += RZ(pi,0)
exec(c)
@test c.density_matrix ≈ [1 0; 0 0]

c = Curcuit(1)
c += H(0)
c += RZ(pi/2,0)
exec(c)
@test c.density_matrix ≈ [.5 -.5im; .5im .5]

# Testing SWAP.

c = Curcuit(2)
c += X(0)
c += SWAP(0,1)
exec(c)
@test c.density_matrix == [0; 0; 1; 0] * [0 0 1 0]

c = Curcuit(4)
c += X(1)
c += X(0)
c += X(3)
c += CCSWAP(0,3,1,2)
exec(c)
@test c.density_matrix == [0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1] * ([0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1])'

c = Curcuit(4)
c += X(1)
c += X(3)
c += CCSWAP(0,3,1,2) # Control-SWAP does nothing because one of the control qubits is 0.
exec(c)
@test c.density_matrix == [0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0] * ([0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0])'

c = Curcuit(2)
c += X(0)
c += H(0)
c += SWAP(0,1)
exec(c)
@test c.density_matrix == [1/√2; 0; -1/√2; 0] * [1/√2 0 -1/√2 0]

c = Curcuit(2)
c += X(0)
c += H(0)
c += SWAP(0,1)
exec(c)
@test c.density_matrix == [1/√2; 0; -1/√2; 0] * [1/√2 0 -1/√2 0]

c = Curcuit(3)
c += X(0)
c += CSWAP(2,0,1)
exec(c)
@test c.density_matrix == [0; 1; 0; 0; 0; 0; 0; 0] * [0 1 0 0 0 0 0 0]

# Testing amplitude damping. Amplitude damping causes |1> to decay to |0> with some probability.

# When we start in |0> state, damping kraus with any |1> to |0> decay probability never gets applied,
# residual kraus always gets applied, so the state doesn't change.
c = Curcuit(1)
c += noisify(ID(0), Gates.damping_kraus_map(.1))
exec(c)
@test c.density_matrix == [1 0; 0 0]

# Starting in |1> state, setting |1> to |0> decay probability to 0.
# Residual kraus should be applied and the state should remain |1>.
c = Curcuit(1)
c += damp_amplitude(X(0), 0.0)
exec(c)
@test c.density_matrix == [0 0; 0 1]

# Starting in |1> state, setting |1> to |0> decay probability to 1.
# Damping kraus should be applied and the state should decay to |0>.
c = Curcuit(1)
c += damp_amplitude(X(0), 1.0)
exec(c)
@test c.density_matrix == [1 0; 0 0]

# Starting in |+> state, setting |1> to |0> decay probability to 0.
# Residual kraus should be applied and the state should remain |1>.
c = Curcuit(1)
c += damp_amplitude(H(0), 0.0)
exec(c)
@test c.density_matrix ≈ [.5 .5; .5 .5]

# Starting in |+> state, setting |1> to |0> decay probability to 1.
# Damping kraus should be applied and the state should decay to |0>.
c = Curcuit(1)
c += damp_amplitude(H(0), 1.0)
exec(c)
@test c.density_matrix ≈ [1 0; 0 0]

# Testing adding control bits to a gate.

p1 = Curcuit(2)
p1 += CNOT(1,0)
exec(p1)
p2 = Curcuit(2)
p2 += controlled(1, X(0))
exec(p2)
@test p1.density_matrix == p2.density_matrix

p1 = Curcuit(3)
p1 += CCNOT(1,2,0)
exec(p1)
p2 = Curcuit(3)
p2 += controlled(1, controlled(2,X(0)))
exec(p2)
@test p1.density_matrix == p2.density_matrix

p1 = Curcuit(3)
p1 += CCNOT(1,2,0)
exec(p1)
p2 = Curcuit(3)
p2 += controlled(1, controlled(2,X(0)))
exec(p2)
@test p1.density_matrix == p2.density_matrix

p1 = Curcuit(3)
p1 += CCNOT(1,2,0)
exec(p1)
p2 = Curcuit(3)
p2 += controlled([1,2], X(0))
exec(p2)
@test p1.density_matrix == p2.density_matrix

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

@test expandGateToFullRegister(1, matrix_H, 0) == matrix_H
@test expandGateToFullRegister(2, HH, 0) == HH
@test expandGateToFullRegister(3, HH, 0, [2]) == CHH
@test expandGateToFullRegister(6, HH, 2, [4]) == matrix_I ⊗ CHH ⊗ matrix_I ⊗ matrix_I
@test expandGateToFullRegister(6, matrix_H, 3, [4,1]) == matrix_I ⊗ CHIC ⊗ matrix_I
@test expandGateToFullRegister(2, matrix_X, 0, [1]) == CX
@test expandGateToFullRegister(3, matrix_X, 0, [1,2]) == CCX
@test expandGateToFullRegister(4, matrix_X, 0, [2,1,3]) == CCCX
@test expandGateToFullRegister(3, matrix_X, 0, [1]) == matrix_I ⊗ CX
@test expandGateToFullRegister(3, matrix_X, 1, [2]) == CX ⊗ matrix_I
@test expandGateToFullRegister(4, matrix_X, 0, [1]) == matrix_I ⊗ matrix_I ⊗ CX
@test expandGateToFullRegister(5, matrix_X, 0, [1,2]) == matrix_I ⊗ matrix_I ⊗ CCX
@test expandGateToFullRegister(2, matrix_X, 1, [0]) == XC
@test expandGateToFullRegister(3, matrix_X, 1, [0]) == matrix_I ⊗ XC
@test expandGateToFullRegister(3, matrix_X, 2, [0,1]) == XCC
@test expandGateToFullRegister(4, matrix_X, 2, [1,0]) == matrix_I ⊗ XCC
@test expandGateToFullRegister(3, matrix_X, 1, [0,2]) == CXC
@test expandGateToFullRegister(4, matrix_X, 2, [1,3]) == CXC ⊗ matrix_I
@test expandGateToFullRegister(5, matrix_X, 2, [3,1]) == matrix_I ⊗ CXC ⊗ matrix_I
@test expandGateToFullRegister(4, matrix_X, 1, [0,3]) == CIXC
@test expandGateToFullRegister(3, matrix_X, 2, [1]) == XC ⊗ matrix_I
@test expandGateToFullRegister(4, matrix_X, 2, [3,0]) == CXIC
@test expandGateToFullRegister(7, matrix_X, 5, [3,6]) == CXIC ⊗ matrix_I ⊗ matrix_I ⊗ matrix_I
@test expandGateToFullRegister(6, matrix_X, 3, [1,4]) == matrix_I ⊗ CXIC ⊗ matrix_I
@test expandGateToFullRegister(1, matrix_S, 0) == matrix_S
#Testing asymmetric gates to ensure the gate doesn't get transposed by mistake.
@test expandGateToFullRegister(1, Gates.damping_residual_kraus(), 0) == Gates.damping_residual_kraus()
@test expandGateToFullRegister(3, Gates.damping_residual_kraus(), 1) == matrix_I ⊗ Gates.damping_residual_kraus() ⊗ matrix_I
@test expandGateToFullRegister(4, Gates.damping_residual_kraus(), 2) == matrix_I ⊗ Gates.damping_residual_kraus() ⊗ matrix_I ⊗ matrix_I
@test expandGateToFullRegister(1, Gates.damping_kraus(), 0) == Gates.damping_kraus()
@test expandGateToFullRegister(3, Gates.damping_kraus(), 1) == matrix_I ⊗ Gates.damping_kraus() ⊗ matrix_I
@test expandGateToFullRegister(4, Gates.damping_kraus(), 1) == matrix_I ⊗ matrix_I ⊗ Gates.damping_kraus() ⊗ matrix_I

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

@test Gates.damping_kraus_map(.1) == [[1 0; 0 √.9], [0 √.1; 0 0]]

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
