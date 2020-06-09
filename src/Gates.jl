module Gates

export Gate, ID, X, Y, Z, H, PHASE, S, T, CZ, RX, RY, RZ, CNOT, CCNOT, SWAP, CSWAP, CCSWAP, custom_gate, controlled, noisify, damp_amplitude

const empty_matrix = Matrix{Complex{Float64}}(undef, 0, 0)
const empty_kraus = Vector{Matrix{Complex{Float64}}}(undef, 0)

struct Gate
	# Gate has either matrix or SWAP bit indexes set. Control bits can be set in either case.
	matrix::Matrix{Complex{Float64}}
	qubit_index::Integer
	control_bit_indexes::AbstractVector{Int64}
	swap_bit_index_a::Integer
	swap_bit_index_b::Integer 
	kraus_operators::AbstractVector{Matrix{Complex{Float64}}}

	# Constructor for all ideal gates except for SWAP.
	function Gate(matrix::Matrix{Complex{Float64}},
		qubit_index::Integer,
		control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))

		return new(matrix, qubit_index, control_bit_indexes, -1, -1, empty_kraus)
	end

	# Constructor for SWAP gates, including controlled once.
	function Gate(swap_bit_index_a::Integer,
		swap_bit_index_b::Integer,
		control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))

		return new(empty_matrix, -1, control_bit_indexes, swap_bit_index_a, swap_bit_index_b, empty_kraus)
	end

	# Constructor for noisy gates.
	function Gate(matrix::Matrix{Complex{Float64}},
		qubit_index::Integer,
		kraus_operators::AbstractVector{Matrix{Complex{Float64}}},
		control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))

		return new(matrix, qubit_index, control_bit_indexes, -1, -1, kraus_operators)
	end

	# Constructor that populates all properties; used to create a new gate based on existing one.
	function Gate(matrix::Matrix{Complex{Float64}},
		qubit_index::Integer,
		control_bit_indexes::AbstractArray{Int64,1},
		swap_bit_index_a::Integer,
		swap_bit_index_b::Integer,
		kraus_operators::AbstractVector{Matrix{Complex{Float64}}})

		return new(matrix, qubit_index, control_bit_indexes, swap_bit_index_a, swap_bit_index_b, kraus_operators)
	end
end

# Identity gate.
# Named ID instead of I to avoid name conflict with LinearAlgebra.I.
function ID(qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([1 0; 0 1]), qubit_index)
end

function X(qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([0 1; 1 0]), qubit_index)
end

function Y(qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([0 complex(0, -1); complex(0, 1) 0]), qubit_index)
end

function Z(qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([1 0; 0 -1]), qubit_index)
end

function H(qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}(1 / √2 * [1 1; 1 -1]), qubit_index)
end

function PHASE(angle::Real, qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([1 0; 0 exp(complex(0, 1) * angle)]), qubit_index)
end

function S(qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([1 0; 0 complex(0, 1)]), qubit_index)
end

function T(qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([1 0; 0 1 / √2 + 1im / √2]), qubit_index)
end

function CZ(control_index::Integer, target_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([1 0; 0 -1]), target_index, [control_index])
end

function RX(angle::Real, qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([cos(angle / 2) complex(0, -1) * sin(angle / 2); complex(0, -1) * sin(angle / 2) cos(angle / 2)]), qubit_index)
end

function RY(angle::Real, qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([cos(angle / 2) -sin(angle / 2); sin(angle / 2) cos(angle / 2)]), qubit_index)
end

function RZ(angle::Real, qubit_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([exp(complex(0, -1) * angle / 2) 0; 0 exp(complex(0, 1) * angle / 2)]), qubit_index)
end

function CNOT(control_index::Integer, target_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([0 1; 1 0]), target_index, [control_index])
end

function CCNOT(control_index_a::Integer, control_index_b::Integer, target_index::Integer)::Gate
	return Gate(Matrix{Complex{Float64}}([0 1; 1 0]), target_index, [control_index_a, control_index_b])
end

function SWAP(swap_bit_index_a::Integer, swap_bit_index_b::Integer)::Gate
	return Gate(swap_bit_index_a, swap_bit_index_b)
end

function CSWAP(control_index::Integer, swap_bit_index_a::Integer, swap_bit_index_b::Integer)::Gate
	return Gate(swap_bit_index_a, swap_bit_index_b, [control_index])
end

function CCSWAP(control_index_a::Integer, control_index_b::Integer, swap_bit_index_a::Integer, swap_bit_index_b::Integer)::Gate
	return Gate(swap_bit_index_a, swap_bit_index_b, [control_index_a, control_index_b])
end

function custom_gate(matrix::Matrix{Complex{Float64}}, qubit_index::Integer, control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([]))::Gate
	return Gate(matrix, qubit_index, control_bit_indexes)
end

function controlled(control_bit_index::Integer, g::Gate)::Gate
	return controlled([control_bit_index], g)
end

function controlled(control_bit_indexes::AbstractArray{Int64,1}, g::Gate)::Gate
	return Gate(g.matrix,
		g.qubit_index,
		vcat(g.control_bit_indexes, control_bit_indexes),
		g.swap_bit_index_a,
		g.swap_bit_index_b,
		g.kraus_operators)
end

function noisify(g::Gate, kraus_operators::AbstractVector{Matrix{Complex{Float64}}})::Gate
	if size(g.kraus_operators)[1] > 0
		error("cannot apply Kraus operators to a gate that already has Kraus operators applied to")
	end
	return Gate(g.matrix,
		g.qubit_index,
		g.control_bit_indexes,
		g.swap_bit_index_a,
		g.swap_bit_index_b,
		kraus_operators)
end

damping_residual_kraus(decay_1_to_0_probability = .1) = Matrix{Complex{Float64}}([1 0; 0 √(1 - decay_1_to_0_probability)])
damping_kraus(decay_1_to_0_probability = .1) = Matrix{Complex{Float64}}([0 √decay_1_to_0_probability; 0 0])
damping_kraus_map(decay_1_to_0_probability = .1) = [damping_residual_kraus(decay_1_to_0_probability), damping_kraus(decay_1_to_0_probability)]

# damp_amplitude adds noise via amplitude damping kraus operators with given probability of decay |1> to |0>.
function damp_amplitude(g::Gate, ket1_to_ket0_decay_probability::Real)
	return noisify(g, damping_kraus_map(ket1_to_ket0_decay_probability))
end

end # module
