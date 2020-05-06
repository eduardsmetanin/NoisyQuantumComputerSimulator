# module NoisyGates

# export noise, amplitude_damping

# include("Gates.jl")

# using .Gates

# # function X(qubit_index::Integer)::Gate
# # 	return Gate(Matrix{Complex{Float64}}([0 1; 1 0]), qubit_index)
# # end

# # noise runs the specified via Kraus operators noise gate.
# function noise(program::Program, 
# 	kraus_operators::AbstractVector{Matrix{Complex{Float64}}},
# 	gate_lowest_index::Integer)::Gate

# 	# println("------------------------------------------------------------------------------------")
# 	# println("kraus_operators: $kraus_operators")
# 	push!(program.commands, ()->begin
# 		kraus_probabilities = Vector{Float64}([])
# 		big_kraus_operators = Vector{Matrix{Complex{Float64}}}([])
# 		for kraus ∈ kraus_operators
# 			big_kraus = FullRegisterGate.build(program.qubit_count, kraus, gate_lowest_index)
# 			# println("big_kraus: $big_kraus")
# 			trace = tr(big_kraus * program.density_matrix * big_kraus')
# 			if imag(trace) ≉ 0
# 				error("expected trace be a real value since it's probability, got $trace")
# 			end
# 			kraus_probability = real(tr(big_kraus * program.density_matrix * big_kraus'))
# 			# println("kraus_probability: $kraus_probability")
# 			push!(big_kraus_operators, big_kraus)
# 			push!(kraus_probabilities, kraus_probability)
# 		end
# 		random_kraus_index = random_index(kraus_probabilities)
# 		# println("random_kraus_index: $random_kraus_index")
# 		random_big_kraus = big_kraus_operators[random_kraus_index]
# 		# println("random_big_kraus: $random_big_kraus")
# 		random_big_kraus_probability = kraus_probabilities[random_kraus_index]
# 		# println("random_big_kraus_probability: $random_big_kraus_probability")
# 		program.density_matrix = random_big_kraus * program.density_matrix * random_big_kraus' / random_big_kraus_probability
# 		# println("program.density_matrix: $(program.density_matrix)")
# 		return
# 	end)
# 	return
# end

# damping_residual_kraus(decay_1_to_0_probability = .1) = Matrix{Complex{Float64}}([1 0; 0 √(1-decay_1_to_0_probability)])
# damping_kraus(decay_1_to_0_probability = .1) = Matrix{Complex{Float64}}([0 √decay_1_to_0_probability; 0 0])
# damping_kraus_map(decay_1_to_0_probability = .1) = [damping_residual_kraus(decay_1_to_0_probability), damping_kraus(decay_1_to_0_probability)]

# # amplitude_damping adds noise via amplitude damping kraus operators with given probability of decay |1> to |0>.
# function amplitude_damping(program::Program, 
# 	ket1_to_ket0_decay_probability::Real,
# 	gate_lowest_index::Integer)

# 	noise(program, damping_kraus_map(ket1_to_ket0_decay_probability), gate_lowest_index)
# 	return
# end

# # random_index returns random index according to the given probabilities.
# # Min returned value is 1; max return value is the length of the given array of probabilities.
# function random_index(probabilities::AbstractVector{Float64})
# 	random = rand()
# 	# println("random: $random")
# 	return value_index_by_probabilities(random, probabilities)
# end

# # value_index_by_probabilities returns index of the given value according to the given probabilities.
# # Min returned value is 1; max return value is the length of the given array of probabilities.
# function value_index_by_probabilities(value::Float64, probabilities::AbstractVector{Float64})
# 	total_probability = zero(Float64)
# 	for (index, probability) ∈ enumerate(probabilities)
# 		total_probability += probability
# 		if value < total_probability
# 			return index
# 		end
# 	end
# 	return length(probabilities) # Unlikely, but control might reach here in case of some rounding issues.
# end

# end # module
