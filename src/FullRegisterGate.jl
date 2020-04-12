module FullRegisterGate

# build expands the given gate with optional control qubits to entire quantum register with the given size.
function build(register_size::Integer,
               gate::AbstractMatrix{Complex{Float64}},
               gate_lowest_index::Integer,
               control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([])
              )::AbstractMatrix{Complex{Float64}}
	small_gate_size = size(gate)[1]
	small_gate_qubit_count = Int(log2(small_gate_size))
	small_gate_bitmask = gate_bitmask(small_gate_qubit_count, gate_lowest_index)
	if small_gate_bitmask == 0
		error("small_gate_bitmask cannot be 0")
	end
	small_gate_offset = first_set_bit_position(small_gate_bitmask)
	if small_gate_offset == -1
		error("internal error: small_gate_offset should not be -1")
	end
	n = 2 ^ register_size
	big_gate = zeros(Complex{Float64}, n, n)
	ctrl_bitmask = control_bitmask(control_bit_indexes)
	for big_gate_row_index ∈ 0:n-1
		if !all_control_bits_set(big_gate_row_index, ctrl_bitmask)
			big_gate[big_gate_row_index+1, big_gate_row_index+1] = Complex(1)
			continue
		end
		target_bits = (big_gate_row_index & small_gate_bitmask) >>> small_gate_offset
		output_state = gate[:,target_bits+1]
		for state_index ∈ 0:small_gate_size-1
			big_gate_column_index = (big_gate_row_index & ~small_gate_bitmask) | (state_index << small_gate_offset)
			big_gate[big_gate_row_index+1, big_gate_column_index+1] = Complex(output_state[state_index+1])
		end
	end
	return big_gate
end

function gate_bitmask(gate_qubit_count::Integer, gate_lowest_index::Integer)::UInt64
	bitmask = zero(UInt64)
	for _ ∈ 1:gate_qubit_count
		bitmask <<= 1
		bitmask |= 1
	end
	return bitmask << gate_lowest_index
end

function control_bitmask(control_bit_indexes::AbstractArray{Int64,1})::UInt64
	bitmask = zero(UInt64)
	for i ∈ control_bit_indexes
		bitmask |= one(UInt64) << i
	end
	return bitmask
end

@inline function all_control_bits_set(i::Integer, control_bitmask::Integer)
	return i & control_bitmask == control_bitmask
end

# first_set_bit_position returns position of the first 1 in the binary
# representation of the given integer counting from the right.
# Returns -1 if 1 not found.
function first_set_bit_position(i::Integer)::Integer
	total_bits = sizeof(i) * 8
	position = 0
	while true
		if i & 1 == 1
			return position
		end
		if position == total_bits - 1
			return -1
		end
		position += 1
		i >>>= 1
	end
end

end # module
