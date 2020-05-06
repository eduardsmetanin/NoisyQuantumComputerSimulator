module FullRegisterGate

export expandGateToFullRegister

# expandGateToFullRegister expands the given gate with optional control qubits to entire quantum register with the given size.
function expandGateToFullRegister(register_size::Integer,
	gate::AbstractMatrix{Complex{Float64}},
	gate_lowest_index::Integer,
	control_bit_indexes::AbstractArray{Int64,1} = Array{Int64,1}([])
	)::AbstractMatrix{Complex{Float64}}

	small_gate_size = size(gate)[1]
	small_gate_qubit_count = Int(log2(small_gate_size))
	gate_bitmask = gate_full_register_bitmask(small_gate_qubit_count, gate_lowest_index)
	if gate_bitmask == 0
		error("gate_bitmask cannot be 0")
	end
	n = 2 ^ register_size
	big_gate = zeros(Complex{Float64}, n, n)
	control_bitmask = control_full_register_bitmask(control_bit_indexes)
	for big_gate_column_index ∈ 0:n-1
		if !all_control_bits_set(big_gate_column_index, control_bitmask)
			big_gate[big_gate_column_index+1, big_gate_column_index+1] = Complex(1)
			continue
		end
		target_bits = (big_gate_column_index & gate_bitmask) >>> gate_lowest_index
		output_state = gate[:,target_bits+1] # Selecting a column here yields the result of matrix-vector multiplication.
		for state_index ∈ 0:small_gate_size-1
			big_gate_row_index = (big_gate_column_index & ~gate_bitmask) | (state_index << gate_lowest_index)
			big_gate[big_gate_row_index+1, big_gate_column_index+1] = Complex(output_state[state_index+1])
		end
	end
	return big_gate
end

function gate_full_register_bitmask(gate_qubit_count::Integer, gate_lowest_index::Integer)::UInt64
	bitmask = zero(UInt64)
	for _ ∈ 1:gate_qubit_count
		bitmask <<= 1
		bitmask |= 1
	end
	return bitmask << gate_lowest_index
end

function control_full_register_bitmask(control_bit_indexes::AbstractArray{Int64,1})::UInt64
	bitmask = zero(UInt64)
	for i ∈ control_bit_indexes
		bitmask |= one(UInt64) << i
	end
	return bitmask
end

@inline function all_control_bits_set(i::Integer, control_bitmask::Integer)
	return i & control_bitmask == control_bitmask
end

end # module
