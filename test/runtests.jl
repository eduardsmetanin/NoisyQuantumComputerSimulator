using Test, NoisyQuantumComputerSimulator, NoisyQuantumComputerSimulator.FullRegisterGate

⊗ = kron

I = Matrix{Complex{Float64}}([1 0; 0 1])
H = Matrix{Complex{Float64}}(1/√2 * [1 1; 1 -1])
X = Matrix{Complex{Float64}}([0 1; 1 0])

# Testing QuantumRegister.

r = QuantumRegister(1)
@test r.density_matrix == [1 0; 0 0]
r = QuantumRegister(2)
@test r.density_matrix == [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
r = QuantumRegister(3)
@test r.density_matrix == [1 0 0 0 0 0 0 0;
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

r = QuantumRegister(2)
h(r,1)
x(r,0,[1])
@test r.density_matrix == BELL₀ * BELL₀'

r = QuantumRegister(2)
x(r,0)
h(r,1)
x(r,0,[1])
@test r.density_matrix == BELL₁ * BELL₁'

r = QuantumRegister(2)
x(r,1)
h(r,1)
x(r,0,[1])
@test r.density_matrix == BELL₂ * BELL₂'

r = QuantumRegister(2)
x(r,0)
x(r,1)
h(r,1)
x(r,0,[1])
@test r.density_matrix == BELL₃ * BELL₃'

r = QuantumRegister(4)
h(r,2)
x(r,1,[2])
@test r.density_matrix ≈ ([1; 0] ⊗ BELL₀ ⊗ [1; 0]) * ([1 0] ⊗ BELL₀' ⊗ [1 0])

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

r = QuantumRegister(3)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₀ * CCBELL₀'

r = QuantumRegister(3)
x(r,0)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₁ * CCBELL₁'

r = QuantumRegister(3)
x(r,1)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₂ * CCBELL₂'

r = QuantumRegister(3)
x(r,0)
x(r,1)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₃ * CCBELL₃'

r = QuantumRegister(3)
x(r,2)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₄ * CCBELL₄'

r = QuantumRegister(3)
x(r,0)
x(r,2)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₅ * CCBELL₅'

r = QuantumRegister(3)
x(r,1)
x(r,2)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₆ * CCBELL₆'

r = QuantumRegister(3)
x(r,0)
x(r,1)
x(r,2)
h(r,0)
h(r,2)
x(r,1,[0,2])
@test r.density_matrix ≈ CCBELL₇ * CCBELL₇'

# Testing h.

r = QuantumRegister(1)
h(r,0)
@test r.density_matrix ≈ [.5 .5; .5 .5]

# Testing x.

r = QuantumRegister(1)
x(r,0)
@test r.density_matrix == [0 0; 0 1]

r = QuantumRegister(1)
gate(r,X,0)
x(r,0)
@test r.density_matrix == [1 0; 0 0]

r = QuantumRegister(1)
h(r,0)
x(r,0)
@test r.density_matrix ≈ [.5 .5; .5 .5]

r = QuantumRegister(1)
gate(r,X,0)
h(r,0)
x(r,0)
@test r.density_matrix ≈ [.5 -.5; -.5 .5]

# Testing y.

r = QuantumRegister(1)
y(r,0)
@test r.density_matrix == [0 0; 0 1]

r = QuantumRegister(1)
x(r,0)
y(r,0)
@test r.density_matrix == [1 0; 0 0]

# Testing z.

r = QuantumRegister(1)
z(r,0)
@test r.density_matrix == [1 0; 0 0]

r = QuantumRegister(1)
x(r,0)
z(r,0)
@test r.density_matrix == [0 0; 0 1]

r = QuantumRegister(1)
h(r,0)
z(r,0)
@test r.density_matrix ≈ [.5 -.5; -.5 .5]

r = QuantumRegister(1)
x(r,0)
h(r,0)
z(r,0)
@test r.density_matrix ≈ [.5 .5; .5 .5]

# Testing swap.

r = QuantumRegister(2)
x(r,0)
swap(r,0,1)
@test r.density_matrix == [0; 0; 1; 0] * [0 0 1 0]

r = QuantumRegister(4)
x(r,1)
x(r,0)
x(r,3)
swap(r,1,2,[0,3])
@test r.density_matrix == [0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1] * ([0; 1] ⊗ [0; 0; 1; 0] ⊗ [0; 1])'

r = QuantumRegister(4)
x(r,1)
x(r,3)
swap(r,1,2,[0,3]) # Control-SWAP does nothing because one of the control qubits is 0.
@test r.density_matrix == [0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0] * ([0; 1] ⊗ [1; 0] ⊗ [0; 1] ⊗ [1; 0])'

r = QuantumRegister(2)
x(r,0)
h(r,0)
swap(r,0,1)
@test r.density_matrix == [1/√2; 0; -1/√2; 0] * [1/√2 0 -1/√2 0]

r = QuantumRegister(2)
x(r,0)
h(r,0)
swap(r,0,1)
@test r.density_matrix == [1/√2; 0; -1/√2; 0] * [1/√2 0 -1/√2 0]

r = QuantumRegister(3)
x(r,0)
swap(r,0,1,[2])
@test r.density_matrix == [0; 1; 0; 0; 0; 0; 0; 0] * [0 1 0 0 0 0 0 0]

# Testing s.

r = QuantumRegister(1)
h(r,0)
s(r,0)
@test r.density_matrix ≈ [.5 -.5im; .5im .5]

r = QuantumRegister(1)
x(r,0)
h(r,0)
s(r,0)
@test r.density_matrix ≈ [.5 .5im; -.5im .5]

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

@test FullRegisterGate.gate_bitmask(1, 0) == 0b1
@test FullRegisterGate.gate_bitmask(1, 1) == 0b10
@test FullRegisterGate.gate_bitmask(1, 2) == 0b100
@test FullRegisterGate.gate_bitmask(2, 0) == 0b11
@test FullRegisterGate.gate_bitmask(2, 1) == 0b110
@test FullRegisterGate.gate_bitmask(2, 2) == 0b1100
@test FullRegisterGate.gate_bitmask(3, 0) == 0b111
@test FullRegisterGate.gate_bitmask(3, 1) == 0b1110
@test FullRegisterGate.gate_bitmask(3, 2) == 0b11100

@test FullRegisterGate.control_bitmask(Array{Int64,1}([])) == 0
@test FullRegisterGate.control_bitmask([0]) == 1
@test FullRegisterGate.control_bitmask([0, 1]) == 0b11
@test FullRegisterGate.control_bitmask([1, 5]) == 0b100010
@test FullRegisterGate.control_bitmask([1, 7, 5]) == 0b10100010

@test FullRegisterGate.all_control_bits_set(0b11, 0b10) == true
@test FullRegisterGate.all_control_bits_set(0b1100, 0b1010) == false
@test FullRegisterGate.all_control_bits_set(0b0000, 0b1010) == false
@test FullRegisterGate.all_control_bits_set(0b1111, 0b0101) == true
@test FullRegisterGate.all_control_bits_set(0b0, 0b0) == true
@test FullRegisterGate.all_control_bits_set(0, 0) == true
@test FullRegisterGate.all_control_bits_set(0xffffffffffffffff, 0xffffffffffffffff) == true
@test FullRegisterGate.all_control_bits_set(0x8fffffffffffffff, 0xffffffffffffffff) == false

@test FullRegisterGate.first_set_bit_position(0b0) == -1
@test FullRegisterGate.first_set_bit_position(0b1) == 0
@test FullRegisterGate.first_set_bit_position(0b10) == 1
@test FullRegisterGate.first_set_bit_position(0b11) == 0
@test FullRegisterGate.first_set_bit_position(0b100) == 2
@test FullRegisterGate.first_set_bit_position(0xffffffffffffffff) == 0
@test FullRegisterGate.first_set_bit_position(0x8000000000000000) == 63
@test FullRegisterGate.first_set_bit_position(0x0000000000000000) == -1

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
