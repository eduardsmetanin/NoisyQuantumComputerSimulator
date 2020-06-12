using Test
using LinearAlgebra
using Curcuits
using Curcuits.FullRegisterGate
using Curcuits.Gates

⊗ = kron

matrix_I = Matrix{Complex{Float64}}([1 0; 0 1])
matrix_H = Matrix{Complex{Float64}}(1 / √2 * [1 1; 1 -1])
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

# Testing reset_state!
c = Curcuit(1, H(0))
@test exec(c) ≈ [.5 .5; .5 .5]
reset_state!(c)
@test c.density_matrix == [1 0; 0 0]

# Testing CNOT.

BELL = 1 / √2 * [1 0 1 0; 0 1 0 1; 0 1 0 -1; 1 0 -1 0] # (CX)(H⊗I)
BELL₀ = BELL * [1; 0; 0; 0]
BELL₁ = BELL * [0; 1; 0; 0]
BELL₂ = BELL * [0; 0; 1; 0]
BELL₃ = BELL * [0; 0; 0; 1]

@test exec(Curcuit(2, H(1), CNOT(1, 0))) == BELL₀ * BELL₀'
@test exec(Curcuit(2, X(0), H(1), CNOT(1, 0))) == BELL₁ * BELL₁'
@test exec(Curcuit(2, X(1), H(1), CNOT(1, 0))) == BELL₂ * BELL₂'
@test exec(Curcuit(2, X(0), X(1), H(1), CNOT(1, 0))) == BELL₃ * BELL₃'
@test exec(Curcuit(4, H(2), CNOT(2, 1))) ≈ ([1; 0] ⊗ BELL₀ ⊗ [1; 0]) * ([1 0] ⊗ BELL₀' ⊗ [1 0])

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

@test exec(Curcuit(3, H(0), H(2), CCNOT(0, 2, 1))) ≈ CCBELL₀ * CCBELL₀'
@test exec(Curcuit(3, X(0), H(0), H(2), CCNOT(0, 2, 1))) ≈ CCBELL₁ * CCBELL₁'
@test exec(Curcuit(3, X(1), H(0), H(2), CCNOT(0, 2, 1))) ≈ CCBELL₂ * CCBELL₂'
@test exec(Curcuit(3, X(0), X(1), H(0), H(2), CCNOT(0, 2, 1))) ≈ CCBELL₃ * CCBELL₃'
@test exec(Curcuit(3, X(2), H(0), H(2), CCNOT(0, 2, 1))) ≈ CCBELL₄ * CCBELL₄'
@test exec(Curcuit(3, X(0), X(2), H(0), H(2), CCNOT(0, 2, 1))) ≈ CCBELL₅ * CCBELL₅'
@test exec(Curcuit(3, X(1), X(2), H(0), H(2), CCNOT(0, 2, 1))) ≈ CCBELL₆ * CCBELL₆'
@test exec(Curcuit(3, X(0), X(1), X(2), H(0), H(2), custom_gate(matrix_X, 1, [0,2]))) ≈ CCBELL₇ * CCBELL₇'

# Testing X.
@test exec(Curcuit(1, X(0))) == [0 0; 0 1]
@test exec(Curcuit(1, H(0), X(0))) ≈ [.5 .5; .5 .5]

# Testing Y.
@test exec(Curcuit(1, Y(0))) == [0 0; 0 1]
@test exec(Curcuit(1, X(0), Y(0))) == [1 0; 0 0]

# Testing Z.
@test exec(Curcuit(1, Z(0))) == [1 0; 0 0]
@test exec(Curcuit(1, X(0), Z(0))) == [0 0; 0 1]
@test exec(Curcuit(1, H(0), Z(0))) ≈ [.5 -.5; -.5 .5]
@test exec(Curcuit(1, X(0), H(0), Z(0))) ≈ [.5 .5; .5 .5]

# Testing H.
@test exec(Curcuit(1, H(0))) ≈ [.5 .5; .5 .5]

# Testing PHASE.
@test exec(Curcuit(1, H(0), PHASE(pi / 2, 0))) ≈ [.5 complex(0, -.5); complex(0, .5) .5]

# Testing S.
@test exec(Curcuit(1, H(0), S(0))) ≈ [.5 -.5im; .5im .5]
@test exec(Curcuit(1, X(0), H(0), S(0))) ≈ [.5 .5im; -.5im .5]

# Testing T.
@test exec(Curcuit(1, H(0), T(0))) ≈ [.5 (1 - 1im) / (2 * sqrt(2)); (1 + 1im) / (2 * sqrt(2)) .5]

# Testing CZ.
@test exec(Curcuit(2, H(0), H(1), CZ(1, 0))) ≈ [.5; .5; .5; -.5] * [.5 .5 .5 -.5]

# Testing RX.
@test exec(Curcuit(1, RX(pi, 0))) ≈ [0 0; 0 1]
@test exec(Curcuit(1, RX(pi / 2, 0))) ≈ [.5 .5im; -.5im .5]

# Testing RY.
@test exec(Curcuit(1, RY(pi, 0))) ≈ [0 0; 0 1]
@test exec(Curcuit(1, RY(pi / 2, 0))) ≈ [.5 .5; .5 .5]

# Testing RZ.
@test exec(Curcuit(1, RZ(pi, 0))) ≈ [1 0; 0 0]
@test exec(Curcuit(1, H(0), RZ(pi, 0))) ≈ [.5 -.5; -.5 .5]
@test exec(Curcuit(1, H(0), RZ(pi / 2, 0))) ≈ [.5 -.5im; .5im .5]

# Testing SWAP.
@test exec(Curcuit(2, X(0), SWAP(0, 1))) == [0; 0; 1; 0] * [0 0 1 0]
@test exec(Curcuit(4, X(1), X(0), X(3), CCSWAP(0, 3, 1, 2))) == [0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1] * ([0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1])'
# Control-SWAP does nothing because one of the control qubits is 0.
@test exec(Curcuit(4, X(1), X(3), CCSWAP(0, 3, 1, 2))) == [0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0] * ([0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0])'
@test exec(Curcuit(2, X(0), H(0), SWAP(0, 1))) == [1 / √2; 0; -1 / √2; 0] * [1 / √2 0 -1 / √2 0]
@test exec(Curcuit(2, X(0), H(0), SWAP(0, 1))) == [1 / √2; 0; -1 / √2; 0] * [1 / √2 0 -1 / √2 0]
@test exec(Curcuit(3, X(0), CSWAP(2, 0, 1))) == [0; 1; 0; 0; 0; 0; 0; 0] * [0 1 0 0 0 0 0 0]

# Testing parametric RX.
c = Curcuit(1, RX("theta", 0))
@test exec(c, Dict("theta" => Float64(pi))) ≈ [0 0; 0 1]
reset_state!(c)
@test exec(c, Dict("theta" => Float64(pi / 2))) ≈ 1 / 2 * [1 1im; -1im 1]

# Testing parametric RY.
c = Curcuit(1, RY("theta", 0))
@test exec(c, Dict("theta" => Float64(pi))) ≈ [0 0; 0 1]
reset_state!(c)
@test exec(c, Dict("theta" => Float64(pi / 2))) ≈ [.5 .5; .5 .5]

# Testing parametric RZ.
c = Curcuit(1, H(0), RZ("theta", 0))
@test exec(c, Dict("theta" => Float64(pi))) ≈ [.5 -.5; -.5 .5]
reset_state!(c)
@test exec(c, Dict("theta" => Float64(pi / 2))) ≈ [.5 -.5im; .5im .5]

# Testing custom_gate without parameters.
@test exec(Curcuit(1, custom_gate(matrix_X, 0), X(0))) == [1 0; 0 0]
@test exec(Curcuit(1, custom_gate(matrix_X, 0), H(0), X(0))) ≈ [.5 -.5; -.5 .5]

# Testing parametric custom_gate.
matrix_x_func = function (params)
    angle = params["alpha"]
    return Matrix{Complex{Float64}}([cos(angle / 2) complex(0, -1) * sin(angle / 2); complex(0, -1) * sin(angle / 2) cos(angle / 2)])
end
c = Curcuit(1, Gate(matrix_x_func, 0))
@test exec(c, Dict("alpha" => Float64(pi))) ≈ [0 0; 0 1]
reset_state!(c)
@test exec(c, Dict("alpha" => Float64(pi / 2))) ≈ 1 / 2 * [1 1im; -1im 1]

# Testing adding control bits to a gate.
@test exec(Curcuit(2, CNOT(1, 0))) == exec(Curcuit(2, controlled(1, X(0))))
@test exec(Curcuit(3, CCNOT(1, 2, 0))) == exec(Curcuit(3, controlled(1, controlled(2, X(0)))))
@test exec(Curcuit(3, CCNOT(1, 2, 0))) == exec(Curcuit(3, controlled(1, controlled(2, X(0)))))
@test exec(Curcuit(3, CCNOT(1, 2, 0))) == exec(Curcuit(3, controlled([1,2], X(0))))

# Testing += sintax for adding gates.
c = Curcuit(2)
c += X(0)
c += X(1)
exec(c)
@test c.density_matrix == exec(Curcuit(2, X(0), X(1)))

# Testing amplitude damping. Amplitude damping causes |1> to decay to |0> with some probability.
# When we start in |0> state, damping kraus with any |1> to |0> decay probability never gets applied,
# residual kraus always gets applied, so the state doesn't change.
@test exec(Curcuit(1, noisify(ID(0), Gates.damping_kraus_map(.1)))) == [1 0; 0 0]
# Starting in |1> state, setting |1> to |0> decay probability to 0.
# Residual kraus should be applied and the state should remain |1>.
@test exec(Curcuit(1, damp_amplitude(X(0), 0.0))) == [0 0; 0 1]
# Starting in |1> state, setting |1> to |0> decay probability to 1.
# Damping kraus should be applied and the state should decay to |0>.
@test exec(Curcuit(1, damp_amplitude(X(0), 1.0))) == [1 0; 0 0]

# Starting in |+> state, setting |1> to |0> decay probability to 0.
# Residual kraus should be applied and the state should remain |1>.
@test exec(Curcuit(1, damp_amplitude(H(0), 0.0))) ≈ [.5 .5; .5 .5]

# Starting in |+> state, setting |1> to |0> decay probability to 1.
# Damping kraus should be applied and the state should decay to |0>.
@test exec(Curcuit(1, damp_amplitude(H(0), 1.0))) ≈ [1 0; 0 0]

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

τ = 1 / √2
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
# Testing asymmetric gates to ensure the gate doesn't get transposed by mistake.
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

@test Curcuits.value_index_by_probabilities(.0, [.3, .4, .2, .1]) == 1
@test Curcuits.value_index_by_probabilities(.35, [.3, .4, .2, .1]) == 2
@test Curcuits.value_index_by_probabilities(.8, [.3, .4, .2, .1]) == 3
@test Curcuits.value_index_by_probabilities(.95, [.3, .4, .2, .1]) == 4
@test Curcuits.value_index_by_probabilities(.9999, [.3, .4, .2, .1]) == 4
@test Curcuits.value_index_by_probabilities(1.0, [.3, .4, .2, .1]) == 4
