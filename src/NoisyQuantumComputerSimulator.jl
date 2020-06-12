module NoisyQuantumComputerSimulator

export Circuit, exec, reset_state!
export Gate, ID, X, Y, Z, H, PHASE, S, T, CZ, RX, RY, RZ, CNOT, CCNOT, SWAP, CSWAP, CCSWAP, controlled, noisify, damp_amplitude
export damping_residual_kraus, damping_kraus, damping_kraus_map

include("FullRegisterGate.jl")
include("Gates.jl")

using LinearAlgebra
using .FullRegisterGate
using .Gates

mutable struct Circuit
	qubit_count::Int
	density_matrix::Matrix{Complex{Float64}}
	commands::Array{Function}

	# Circuit creates quantum register with the given size and sets it's value to |0⟩.
	function Circuit(size::Integer, gates...)
		state_size = 2^size
		density_matrix = zeros(Complex{Float64}, state_size, state_size)
		density_matrix[1,1] = 1
		Circuit = new(Int(size), density_matrix, Array{Function}[])
		for gate ∈ gates
			Circuit += gate
		end
		return Circuit
	end
end

# reset_state! sets quantum register state to |0⟩.
function reset_state!(c::Circuit)
	c.density_matrix .= complex(0.0, 0.0)
	c.density_matrix[1,1] = 1
end

function add_gate(c::Circuit, g::Gate)::Circuit
	if !isempty(g.control_bit_indexes) && !isempty(g.kraus_operators)
		error("applying Kraus operators within a gate with control bit(s) is not supported; instead, apply Kraus operators to control and target qubits individually")
	end
	for control ∈ g.control_bit_indexes
		if control < 0 || control >= c.qubit_count
			error("control bit at index $control is outside of bounds of quantum register [0, $(c.qubit_count - 1)]")
		end
	end
	if g.swap_bit_index_a ≠ -1 && g.swap_bit_index_b ≠ -1
		if g.swap_bit_index_a < 0 || g.swap_bit_index_a >= c.qubit_count
			error("bit to SWAP at index $(g.swap_bit_index_a) is outside of bounds of quantum register [0, $(c.qubit_count - 1)]")
		end
		if g.swap_bit_index_b < 0 || g.swap_bit_index_b >= c.qubit_count
			error("bit to SWAP at index $(g.swap_bit_index_b) is outside of bounds of quantum register [0, $(c.qubit_count - 1)]")
		end
		# build SWAP as a sequence of 3 CNOTs.
		not_matrix = Matrix{Complex{Float64}}([0 1; 1 0])
		controls_a = vcat(g.swap_bit_index_a, g.control_bit_indexes)
		controls_b = vcat(g.swap_bit_index_b, g.control_bit_indexes)
		add_unitary_gate!(c, not_matrix, ()->(), g.swap_bit_index_a, controls_b)
		add_unitary_gate!(c, not_matrix, ()->(), g.swap_bit_index_b, controls_a)
		add_unitary_gate!(c, not_matrix, ()->(), g.swap_bit_index_a, controls_b)
	else
		add_unitary_gate!(c, g.matrix, g.matrix_func, g.qubit_index, g.control_bit_indexes)
	end
	if size(g.kraus_operators)[1] > 0
		add_noise!(c, g.kraus_operators, g.qubit_index)
	end
	return c
end

# Defining + allows adding gates like this: Circuit += gate
Base.:+(c::Circuit, g::Gate) = add_gate(c, g)

function exec(c::Circuit, params::Dict{String,Float64} = Dict{String,Float64}())::Array{Complex{Float64}}
	for command ∈ c.commands
		command(params)
	end
	return c.density_matrix
end

# add_unitary_gate! adds the given unitary gate with optional control bits to the given Circuit.
function add_unitary_gate!(c::Circuit, 
	matrix::AbstractMatrix{Complex{Float64}},
	matrix_func::Function,
	gate_lowest_index::Integer,
	control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))

	push!(c.commands, (params::Dict{String,Float64})->begin
		small_matrix = !isempty(matrix) ? matrix : matrix_func(params)
		matrix_size = size(small_matrix)[1]
		if matrix_size == 0
			error("matrix is empty")
		end
		gate_size = Int(log2(matrix_size))
		if gate_lowest_index < 0 || gate_lowest_index + gate_size > c.qubit_count
			error("quantum register of size $(c.qubit_count) cannot fit given gate of size $(gate_size) at index $(gate_lowest_index)")
		end
		big_matrix = expandGateToFullRegister(c.qubit_count, small_matrix, gate_lowest_index, control_bit_indexes)
		c.density_matrix = big_matrix * c.density_matrix * big_matrix'
		return
	end)
	return
end

# add_noise! adds the specified via Kraus operators noise gate to the given Circuit.
function add_noise!(c::Circuit, 
	kraus_operators::AbstractVector{Matrix{Complex{Float64}}},
	gate_lowest_index::Integer)

	push!(c.commands, (params::Dict{String,Float64})->begin
		kraus_probabilities = Vector{Float64}([])
		big_kraus_operators = Vector{Matrix{Complex{Float64}}}([])
		for kraus ∈ kraus_operators
			big_kraus = expandGateToFullRegister(c.qubit_count, kraus, gate_lowest_index)
			trace = tr(big_kraus * c.density_matrix * big_kraus')
			if imag(trace) ≉ 0
				error("expected trace be a real value since it's a probability, got $trace")
			end
			kraus_probability = real(tr(big_kraus * c.density_matrix * big_kraus'))
			push!(big_kraus_operators, big_kraus)
			push!(kraus_probabilities, kraus_probability)
		end
		random_kraus_index = random_index(kraus_probabilities)
		big_kraus = big_kraus_operators[random_kraus_index]
		big_kraus_probability = kraus_probabilities[random_kraus_index]
		c.density_matrix = big_kraus * c.density_matrix * big_kraus' / big_kraus_probability
		return
	end)
	return
end

# random_index returns random index according to the given probabilities.
# Min returned value is 1; max return value is the length of the given array of probabilities.
function random_index(probabilities::AbstractVector{Float64})
	random = rand()
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

end # module
