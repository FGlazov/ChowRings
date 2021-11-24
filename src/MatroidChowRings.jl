module MatroidChowRings

using Oscar;
const pm = Polymake;

struct MChowRing
    chow_ring # TODO Add types here.
    projection
    indeterminates::Vector{fmpq_mpoly}
    matroid::pm.BigObject
end

function direct_sum_decomp(matroid::pm.BigObject, matroid_element::Int64)
    direct_sum_decomp(matroid_chow_ring(matroid), matroid_element)
end

function direct_sum_decomp(chow_ring::MChowRing, matroid_element::Int64)
    # TODO: Check bounds of matroid element?

    # TODO: matroid.COLOOPS seems to be broken. Ask why?
    matroid = chow_ring.matroid

    first_term = matroid_chow_ring(pm.matroid.deletion(matroid, matroid_element))

    ground_set = Polymake.Set(range(0, matroid.N_ELEMENTS - 1,step=1))
    flats = matroid.LATTICE_OF_FLATS.FACES
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)
    proper_flats = flats[2:pm.size(flats, 1)-1, 1:pm.size(flats,2)]
    n_proper_flats = pm.size(proper_flats, 1)

    flat_canidates = []
    # TODO: I'm sure this can be done MUCH faster using the IncidenceMatrix stucture.
    # That way, we can turn a "O(2*n)" pass into a "O(n) + O(log n)" pass
    for flat_index in 1:n_proper_flats
        proper_flat = pm.row(proper_flats, flat_index)
        if Polymake.in(matroid_element, proper_flat)
            push!(flat_canidates, proper_flat)
        end
    end
    flat_canidates = Set(flat_canidates)

    # TODO: To make this quicker, we could parallelize this task.
    second_term = []
    for flat_index in 1:n_proper_flats
        proper_flat = pm.row(proper_flats, flat_index)
        if !Polymake.in(matroid_element, proper_flat)
            proper_flat_copy = copy(proper_flat)
            push!(proper_flat_copy, matroid_element)
            if proper_flat_copy in flat_canidates
                matroid_1 = pm.matroid.contraction(matroid, proper_flat_copy)
                matroid_2 = pm.matroid.deletion(matroid, set_complement(ground_set, proper_flat_copy))

                push!(second_term, (matroid_chow_ring(matroid_1), matroid_chow_ring(matroid_2)))
            end
        end
    end

    coloops = matroid.DUAL.LOOPS
    is_coloop = pm.in(matroid_element, coloops)
    coloop_term = nothing
    if is_coloop
        coloop_term = first_term
    end

    #TODO Also return a function which represents the cannonical isomorphism.

    first_term, second_term, coloop_term

end

function matroid_chow_ring(matroid::pm.BigObject)::MChowRing
    if pm.type_name(matroid) != "Matroid"
        throw(ArgumentError("BigObject is not a matroid."))
    end

    flats = matroid.LATTICE_OF_FLATS.FACES;
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)

    # Technically non-empty proper flats. First element is either empty or the whole set, last element is the other of the two.
    proper_flats = flats[2:pm.size(flats, 1)-1, 1:pm.size(flats,2)]

    if length(proper_flats) == 0
        # TODO: Better representation of trivial ring?
        return MChowRing(nothing, nothing, [], matroid)
    end


    base_ring, indeterminates = generate_base_ring(proper_flats)


    type_i_polynomials = generate_type_i_ideal(base_ring, proper_flats, indeterminates, matroid)
    type_j_polynomials = generate_type_j_ideal(proper_flats, indeterminates)
    generators = vcat(type_i_polynomials, type_j_polynomials)

    chow_modulus = ideal(base_ring, generators)
    # Add std(ideal) to compute standard basis?
    chow_ring, projection = quo(base_ring, chow_modulus)

    MChowRing(chow_ring, projection, indeterminates, matroid)
end


function augmented_matroid_chow_ring(matroid::pm.BigObject)
    if pm.type_name(matroid) != "Matroid"
        throw(ArgumentError("BigObject is not a matroid."))
    end

    n_elements = matroid.N_ELEMENTS
    flats = matroid.LATTICE_OF_FLATS.FACES
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)

    # Technically non-empty proper flats. First element is either empty or the whole set, last element is the other of the two.
    proper_flats = flats[2:pm.size(flats, 1)-1, 1:pm.size(flats,2)]
    top_node = matroid.LATTICE_OF_FLATS.TOP_NODE
    if top_node == 0
        proper_flats = flats[2:pm.size(flats, 1), 1:pm.size(flats,2)]
    else
        proper_flats = flats[1:pm.size(flats, 1)-1, 1:pm.size(flats,2)]
    end

    base_ring, matroid_element_vars, flat_vars = generate_augmented_base_ring(proper_flats, n_elements)
    type_i_polynomials = generate_augmented_type_i_ideal(proper_flats, matroid_element_vars, flat_vars)
    type_j_polynomials = generate_augmented_type_j_ideal(proper_flats, matroid_element_vars, flat_vars)
    generators = vcat(type_i_polynomials, type_j_polynomials)

    chow_modulus = ideal(base_ring, generators)
    # Add std(ideal) to compute standard basis?
    chow_ring, projection = quo(base_ring, chow_modulus)
    chow_ring, projection, matroid_element_vars, flat_vars
end


function generate_base_ring(proper_flats)
    variable_names = create_flat_variables_names(proper_flats)
    # TODO: This doesn't seem to work well when it represents the trivial ring.
    #       Printing a trivial ring generated here causes errors.
    #       Perhaps a different way to create the trivial ring?
    PolynomialRing( QQ, variable_names);
end

function generate_augmented_base_ring(proper_flats, n_elements)
    matroid_element_variable_names = Array{String}(undef, n_elements)
    for i in 1:n_elements
        matroid_element_variable_names[i] = "y_" * string(i)
    end

    flat_variable_names = create_flat_variables_names(proper_flats)
    variable_names = vcat(matroid_element_variable_names, flat_variable_names)

    ring, variables = PolynomialRing( QQ, variable_names);

    ring, variables[1:n_elements], variables[n_elements + 1: length(variables)]

end

function create_flat_variables_names(flats)
    n_flats = pm.size(flats, 1)
    variable_names = Array{String}(undef, n_flats)

    for i in 1:n_flats
        row = pm.row(flats, i)
        # TODO: Is there some better way to extract the polymake string without the type?
        set_name = split(string(pm.row(flats, i)), "\n")[2]
        variable_descriptor = replace(set_name[2:length(set_name)-1], " " => "_")
        variable_names[i] = "x__" * variable_descriptor
    end

    variable_names
end


function generate_type_i_ideal(base_ring, proper_flats, indeterminates, matroid)
    n = matroid.N_ELEMENTS
    n_proper_flats = pm.size(proper_flats, 1)
    ideal_polynomials = Vector{fmpq_mpoly}()

    for i in 1:n
        for j in i+1:n
            ij_polynomial = base_ring() # Create zero.
            ij_set = Set{Int64}([i,j])
            for flat_index in 1:n_proper_flats
                proper_flat = pm.row(proper_flats, flat_index)
                if pm.in(i, proper_flat)
                    ij_polynomial += indeterminates[flat_index]
                end
                if pm.in(j, proper_flat)
                    ij_polynomial -= indeterminates[flat_index]
                end
            end
            push!(ideal_polynomials, ij_polynomial)
        end
    end

    ideal_polynomials
end

function generate_augmented_type_i_ideal(proper_flats, matroid_element_vars, flat_vars)
    n_proper_flats = pm.size(proper_flats, 1)
    ideal_polynomials = Array{fmpq_mpoly}(undef, length(matroid_element_vars))

    for i in 1:length(matroid_element_vars)
        ideal_polynomials[i] = matroid_element_vars[i]
        for j in 1:n_proper_flats
            flat = Polymake.row(proper_flats, j)
            # TODO: There must be some much quicker way to access this.
            # I.e. some way to get all rows quickly which don't contain a given element.
            if !Polymake.in(i, flat)
                ideal_polynomials[i] -= flat_vars[j]
            end
        end
    end

    ideal_polynomials
end

function generate_type_j_ideal(proper_flats, indeterminates)
    n_proper_flats = pm.size(proper_flats, 1)
    ideal_polynomials = Vector{fmpq_mpoly}()

    for i in 1:n_proper_flats
        for j in i:n_proper_flats
            i_flat = Polymake.row(proper_flats, i)
            j_flat = Polymake.row(proper_flats, j)
            if are_sets_incomparable(i_flat, j_flat)
                polynomial = indeterminates[i] * indeterminates[j]
                push!(ideal_polynomials, polynomial)
            end
        end
    end

    ideal_polynomials

end


function generate_augmented_type_j_ideal(proper_flats, matroid_element_vars, flat_vars)
    n_proper_flats = pm.size(proper_flats, 1)

    incomparable_polynomials = generate_type_j_ideal(proper_flats, flat_vars)

    xy_polynomials = Vector{fmpq_mpoly}()
    for i in 1:n_proper_flats
        flat = Polymake.row(proper_flats, i)
        for j in 1:length(matroid_element_vars)
            # TODO: There must be some much quicker way to access this.
            # I.e. some way to get all rows quickly which don't contain a given element.
            if !Polymake.in(j, flat)
                push!(xy_polynomials, matroid_element_vars[j] * flat_vars[i])
            end
        end
    end

    vcat(incomparable_polynomials, xy_polynomials)

end

function are_sets_incomparable(a, b)
    return !is_subset(a,b) && !is_subset(b,a)
end

function are_sets_equal(a,b)


    is_subset(a,b) && is_subset(b,a)
end

function is_subset(a, b)
    for element in a
        if !Polymake.in(element, b)
            return false
        end
    end
    true
end

function set_complement(ground_set, to_complement)
    result = copy(ground_set)
    for i in to_complement
        delete!(result, i)
    end
    result
end

end
