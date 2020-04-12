module NoisyQuantumComputerSimulator

export QuantumRegister, gate, h, x, y, z, swap, s

include("FullRegisterGate.jl")

mutable struct QuantumRegister
	qubit_count:: Int
	density_matrix::Array{Complex{Float64}}

	# QuantumRegister creates quantum register with the given size and sets it's value to |0⟩
	function QuantumRegister(size::Integer)
		state_size = 2^size
		density_matrix = zeros(Complex{Float64}, state_size, state_size)
		density_matrix[1,1] = 1
		return new(Int(size), density_matrix)
	end
end

# gate runs the given gate on the given quantum register with optional control bits.
function gate(register::QuantumRegister, 
              matrix::AbstractMatrix{Complex{Float64}},
              gate_lowest_index::Integer,
	      control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]) # TODO: Can it be Integer?
	     )
	big_matrix = FullRegisterGate.build(register.qubit_count, matrix, gate_lowest_index, control_bit_indexes)
	register.density_matrix = big_matrix * register.density_matrix * big_matrix'
	return
end

# gate runs Hadamard gate on the given quantum register with optional control bits.
function h(register::QuantumRegister, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(register, Matrix{Complex{Float64}}(1/√2 * [1 1; 1 -1]), index, control_bit_indexes)
	return
end

# gate runs X gate on the given quantum register with optional control bits.
function x(register::QuantumRegister, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(register, Matrix{Complex{Float64}}([0 1; 1 0]), index, control_bit_indexes)
	return
end

# gate runs Y gate on the given quantum register with optional control bits.
function y(register::QuantumRegister, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(register, Matrix{Complex{Float64}}([0 complex(0, -1); complex(0, 1) 0]), index, control_bit_indexes)
	return
end

# gate runs Z gate on the given quantum register with optional control bits.
function z(register::QuantumRegister, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(register, Matrix{Complex{Float64}}([1 0; 0 -1]), index, control_bit_indexes)
	return
end

# gate runs SWAP gate on the given quantum register with optional control bits.
function swap(register::QuantumRegister,
              index_a::Integer,
              index_b::Integer,
              control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([])
             )
	controls_a = vcat(index_a, control_bit_indexes)
	controls_b = vcat(index_b, control_bit_indexes)
	x(register, index_a, controls_b)
	x(register, index_b, controls_a)
	x(register, index_a, controls_b)
	return
end

# gate runs phase gate (S) on the given quantum register with optional control bits.
function s(register::QuantumRegister, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(register, Matrix{Complex{Float64}}([1 0; 0 complex(0,1)]), index, control_bit_indexes)
	return
end

# function product!(state::AbstractArray, gate::AbstractMatrix{Complex{Float64}})
# 	s = size(state)
# 	if s[1] == 1
# 		error("state cannot be a row vector; use column vector or density matrix")
# 	end
# 	if length(s) == 1 # state is a column (e.g. size([10; 20]) == (2,))
# 		return gate * state
# 	end
# 	return gate * state * gate'
# end

end # module
