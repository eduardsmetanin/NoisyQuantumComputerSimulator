module NoisyQuantumComputerSimulator

using LinearAlgebra

export Program, exec, gate, h, x, y, z, swap, s, noise, amplitude_damping

include("FullRegisterGate.jl")

mutable struct Program
	qubit_count:: Int
	density_matrix::Array{Complex{Float64}}
	commands::Array{Function}

	# Program creates quantum register with the given size and sets it's value to |0⟩.
	function Program(size::Integer)
		state_size = 2^size
		density_matrix = zeros(Complex{Float64}, state_size, state_size)
		density_matrix[1,1] = 1
		return new(Int(size), density_matrix, Array{Function}[])
	end
end

function exec(program::Program)
	for command ∈ program.commands
		command()
	end
	return
end

# gate runs the given gate with optional control bits.
function gate(program::Program, 
	matrix::AbstractMatrix{Complex{Float64}},
	gate_lowest_index::Integer,
	control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]) # TODO: Can it be Integer?
	)
	push!(program.commands, ()->begin
		big_matrix = FullRegisterGate.build(program.qubit_count, matrix, gate_lowest_index, control_bit_indexes)
		program.density_matrix = big_matrix * program.density_matrix * big_matrix'
		return
	end)
	return
end

# h runs Hadamard gate with optional control bits.
function h(program::Program, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(program, Matrix{Complex{Float64}}(1/√2 * [1 1; 1 -1]), index, control_bit_indexes)
	return
end

# x runs X gate with optional control bits.
function x(program::Program, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(program, Matrix{Complex{Float64}}([0 1; 1 0]), index, control_bit_indexes)
	return
end

# y runs Y gate with optional control bits.
function y(program::Program, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(program, Matrix{Complex{Float64}}([0 complex(0, -1); complex(0, 1) 0]), index, control_bit_indexes)
	return
end

# z runs Z gate with optional control bits.
function z(program::Program, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(program, Matrix{Complex{Float64}}([1 0; 0 -1]), index, control_bit_indexes)
	return
end

# swap runs SWAP gate with optional control bits.
function swap(program::Program,
	index_a::Integer,
	index_b::Integer,
	control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))

	controls_a = vcat(index_a, control_bit_indexes)
	controls_b = vcat(index_b, control_bit_indexes)
	x(program, index_a, controls_b)
	x(program, index_b, controls_a)
	x(program, index_a, controls_b)
	return
end

# s runs phase gate (S) with optional control bits.
function s(program::Program, index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))
	gate(program, Matrix{Complex{Float64}}([1 0; 0 complex(0,1)]), index, control_bit_indexes)
	return
end

# noise runs the specified via Kraus operators noise gate.
function noise(program::Program, 
	kraus_operators::AbstractVector{Matrix{Complex{Float64}}},
	gate_lowest_index::Integer)

	# println("------------------------------------------------------------------------------------")
	# println("kraus_operators: $kraus_operators")
	push!(program.commands, ()->begin
		kraus_probabilities = Vector{Float64}([])
		big_kraus_operators = Vector{Matrix{Complex{Float64}}}([])
		for kraus ∈ kraus_operators
			big_kraus = FullRegisterGate.build(program.qubit_count, kraus, gate_lowest_index)
			# println("big_kraus: $big_kraus")
			trace = tr(big_kraus * program.density_matrix * big_kraus')
			if imag(trace) ≉ 0
				error("expected trace be a real value since it's probability, got $trace")
			end
			kraus_probability = real(tr(big_kraus * program.density_matrix * big_kraus'))
			# println("kraus_probability: $kraus_probability")
			push!(big_kraus_operators, big_kraus)
			push!(kraus_probabilities, kraus_probability)
		end
		random_kraus_index = random_index(kraus_probabilities)
		# println("random_kraus_index: $random_kraus_index")
		random_big_kraus = big_kraus_operators[random_kraus_index]
		# println("random_big_kraus: $random_big_kraus")
		random_big_kraus_probability = kraus_probabilities[random_kraus_index]
		# println("random_big_kraus_probability: $random_big_kraus_probability")
		program.density_matrix = random_big_kraus * program.density_matrix * random_big_kraus' / random_big_kraus_probability
		# println("program.density_matrix: $(program.density_matrix)")
		return
	end)
	return
end

damping_residual_kraus(decay_1_to_0_probability = .1) = Matrix{Complex{Float64}}([1 0; 0 √(1-decay_1_to_0_probability)])
damping_kraus(decay_1_to_0_probability = .1) = Matrix{Complex{Float64}}([0 √decay_1_to_0_probability; 0 0])
damping_kraus_operators(decay_1_to_0_probability = .1) = [damping_residual_kraus(decay_1_to_0_probability), damping_kraus(decay_1_to_0_probability)]

# amplitude_damping adds noise via amplitude damping kraus operators with given probability of decay |1> to |0>.
function amplitude_damping(program::Program, 
	ket1_to_ket0_decay_probability::Real,
	gate_lowest_index::Integer)

	noise(program, damping_kraus_operators(ket1_to_ket0_decay_probability), gate_lowest_index)
	return
end

# random_index returns random index according to the given probabilities.
# Min returned value is 1; max return value is the length of the given array of probabilities.
function random_index(probabilities::AbstractVector{Float64})
	random = rand()
	# println("random: $random")
	return value_index_by_probabilities(random, probabilities)
end

# value_index_by_probabilities returns index of the given value according to the given probabilities.
# Min returned value is 1; max return value is the length of the given array of probabilities.
function value_index_by_probabilities(value::Float64, probabilities::AbstractVector{Float64})
	total_probability = zero(Float64)
	for (index, probability) ∈ enumerate(probabilities)
		total_probability += probability
		if value < total_probability
			return index
		end
	end
	return length(probabilities) # Unlikely, but control might reach here in case of some rounding issues.
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
