MultivariateDiscreteNonParametric(xs, ps; map_func=ReversibleMapper, kwargs...) = DiscreteNonParametric(map_func.(xs), ps; kwargs...)

# Define the ReversibleMapper struct
struct ReversibleMapper
    categories        # Categories for each dimension
    indices::Vector{Dict}         # Mapping from category to index
    strides::Vector{Int}                   # Strides for each dimension

    # Constructor to initialize the mapper
    function ReversibleMapper(categories)
        # Create dictionaries mapping each category to its index
        indices = [Dict(cat => idx for (idx, cat) in enumerate(dim)) for dim in categories]

        # Calculate strides based on the size of each dimension
        strides = Vector{Int}(undef, length(categories))
        strides[1] = 1
        for i in 2:length(categories)
            strides[i] = strides[i-1] * length(categories[i-1])
        end

        new(categories, indices, strides)
    end
end


# Function to encode a tuple into a single index
function encode(mapper::ReversibleMapper, tuple::Vector)
    idx = 0
    for (i, val) in enumerate(tuple)
        if !haskey(mapper.indices[i], val)
            error("Value '$val' not found in dimension $i.")
        end
        # Subtract 1 to make it zero-based
        idx += (mapper.indices[i][val] - 1) * mapper.strides[i]
    end
    return idx
end

# Corrected function to decode an index back into a tuple
function decode(mapper::ReversibleMapper, index::Int)
    tuple = Vector(undef, length(mapper.categories))
    for i in reverse(1:length(mapper.categories))
        stride = mapper.strides[i]
        dim_size = length(mapper.categories[i])
        idx_in_dim = (index ÷ stride) % dim_size
        tuple[i] = mapper.categories[i][idx_in_dim+1]
    end
    return tuple
end

function test_example()
    # Define categories for each dimension
    categories = [
        ["Red", "Green", "Blue"],      # Dimension 1
        ["Small", "Medium", "Large"],  # Dimension 2
        ["Circle", "Square"]            # Dimension 3
    ]

    # Initialize the mapper
    mapper = ReversibleMapper(categories)

    # Define multiple tuples to test
    test_tuples = [
        ["Red", "Small", "Circle"],
        ["Green", "Medium", "Square"],
        ["Blue", "Large", "Circle"],
        ["Blue", "Medium", "Square"],
        ["Red", "Large", "Square"]
    ]

    for original_tuple in test_tuples
        println("\nOriginal Tuple: ", original_tuple)

        # Encode the tuple
        index = encode(mapper, original_tuple)
        println("Encoded Index: ", index)

        # Decode back to tuple
        decoded_tuple = decode(mapper, index)
        println("Decoded Tuple: ", decoded_tuple)

        # Verify reversibility
        @assert original_tuple == decoded_tuple
        println("Mapping is reversible!")
    end

    # Additional Test: Ensure all possible combinations are unique and reversible
    println("\nTesting all possible combinations for uniqueness and reversibility:")
    total_combinations = prod(length.(categories))
    unique_indices = Set{Int}()
    for i in 1:total_combinations
        # Generate tuple from index
        decoded = decode(mapper, i - 1)
        encoded = encode(mapper, decoded)
        @assert (i - 1) == encoded
        push!(unique_indices, encoded)
    end
    println("All ", total_combinations, " combinations are unique and reversible.")
end