module MatroidChowRings

export direct_sum_decomp, matroid_chow_ring, augmented_matroid_chow_ring

using Oscar;
const pm = Polymake;

"""
Represents a matroid Chow ring. The result of a matroid_chow_ring call.

Fields
- chow_ring:      The Chow ring of the matroid, represented as an affine algebra.
                  It may be == None, in that case chow ring is trivial (= 0).
- indeterminates: The indeterminates of the chow ring. These are of the form
                  x_F, where F is a nonempty proper flat. The name of F
                  corresponds to elements inside that flat. E.g. x__1_5_7 would
                  be the flat consisting of the elements 1, 5, and 7. Another
                  way to see all the flats would be to call
                  matroid.LATTICE_OF_FLATS.FACES.
                  This vector may be empty, in that case the chow ring is trivial.
- matroid:        The matroid whose chow ring has been computed. It is
                  represented as a polymake Matroid.
"""
struct MChowRing
    chow_ring # TODO Add types here.
    indeterminates::Vector{MPolyQuoElem{fmpq_mpoly}}
    matroid::pm.BigObject
end

"""
Represents a homorphism which sends the RHS of a Chow ring direct sum
decomposition to the LHS of the decompositon.

Note that the right hand side of the decomposition is not represented as
subalgebras of CH(M), but as Chow rings which are isomorphic to the subalgebras
of CH(M). If you wish to turn an element on the RHS to one which lives in the
domain of the LHS, you need to use this structure.

This structure represents the morphism:

first_term_morphism + coloop_term_morphism + sum_(f in second_term_morphism) f( g_f * h_f)

Where g_f,h_f are projections from CH(M_F+1) and CH(M^F), respectively, to CH(M-i).
The input to g_f and h_f would then be values in the corresponding tuple in MChowRingDecomp.

Fields:
- deleted_element:      The element deleted on the RHS of the decomposition.
- first_term_morphism:  Sends elements of the first_term into the LHS.
- second_term_morphism: A vector of morphisms. Sends an element of
- coloop_term_morphism: An optional morphism, that sends members of the
                        coloopterm into the LHS.
"""
struct MChowRingHomorphism
    deleted_element::Int64
    first_term_morphism
    second_term_morphism
    coloop_term_morphism
end

"""
Represents a direct sum decomposition of the Chow ring of a matroid. See
"A semi-small decomposition of the Chow ring of a matroid by Braden et. al" for
details. The result of a direct_sum_decomp call.

Note that the right hand side of the decomposition is not represented as
subalgebras of CH(M), but as Chow rings which are isomorphic to the subalgebras
of CH(M). If you wish to turn an element on the RHS to one which lives in the
domain of the LHS, you need to use the homomorphism field.

Fields:
- deleted_element: The element which was deleted on the RHS of the decomposition.
- chow_ring_LHS:   The chow ring which was decomposed.
- first_term:      Always present part of the decomposition, equal to the chow
                   ring of M-i.
- second_term:     A vector of tuples of chow rings. This has varying length,
                   and each tuple is of the form (CH(M_(F+i)), CH(M^F), where
                   i is the deleted element, and F is a flat.
- coloopterm:      A term that appears in the decomposition if deleted_element
                   is a coloop. Equal to the chow ring of M-i if present.
- homomorphism:    Represents the homomorphism which takes a combination of
                   elements from the first/second_term/coloop_term and maps it
                   to the LHS of the decomposition.
"""
struct MChowRingDecomp #TODO Add types here.
    deleted_element::Int64
    chow_ring_LHS::MChowRing
    first_term::MChowRing
    second_term
    coloopterm
    homomorphism::MChowRingHomorphism
end

#TODO: Document how to use this.


function direct_sum_decomp(matroid::pm.BigObject, matroid_element::Int64)
    direct_sum_decomp(matroid_chow_ring(matroid), matroid_element)
end

function direct_sum_decomp(chow_ring::MChowRing, matroid_element::Int64)
    # TODO: Check bounds of matroid element.
    matroid = chow_ring.matroid

    first_term = matroid_chow_ring(pm.matroid.deletion(matroid, matroid_element))
    first_term_morphism = create_theta_i(first_term, chow_ring, matroid_element)

    ground_set = Polymake.Set(range(0, matroid.N_ELEMENTS - 1,step=1))
    flats = matroid.LATTICE_OF_FLATS.FACES
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)
    proper_flats = flats[2:pm.size(flats, 1)-1, 1:pm.size(flats,2)]
    n_proper_flats = pm.size(proper_flats, 1)

    flat_canidates = []
    # TODO: I'm sure this can be done MUCH faster using the IncidenceMatrix stucture.
    # That way, we can turn a "O(2*n)" pass into a "O(n) + O(log n)" pass
    for flat_index in 1:n_proper_flats
        proper_flat = extract_flat(proper_flats, flat_index)
        if Polymake.in(matroid_element, proper_flat)
            push!(flat_canidates, proper_flat)
        end
    end
    flat_canidates = Set(flat_canidates)

    # This loop might be a good candiate for parralelization.
    second_term = []
    second_term_morphism = []
    for flat_index in 1:n_proper_flats
        proper_flat = extract_flat(proper_flats, flat_index)
        if !Polymake.in(matroid_element, proper_flat)
            proper_flat_copy = copy(proper_flat)
            push!(proper_flat_copy, matroid_element)
            if proper_flat_copy in flat_canidates
                matroid_1 = pm.matroid.contraction(matroid, proper_flat_copy)
                matroid_2 = pm.matroid.deletion(matroid, set_complement(ground_set, proper_flat))

                chow_1 = matroid_chow_ring(matroid_1)
                chow_2 = matroid_chow_ring(matroid_2)

                push!(second_term, (chow_1, chow_2))

                # TODO: Need to adjust term for reindexing here
                # term = find_flat_variable(first_term, proper_flat_copy)
                # morphism = create_theta_i(coloop_term, chow_ring, matroid_element, term)
            end
        end
    end

    coloops = pm.matroid.coloops(matroid)
    is_coloop = pm.in(matroid_element, coloops)
    coloop_term = nothing
    coloop_term_morphism = nothing
    if is_coloop
        coloop_term = first_term
        flat = copy(ground_set)
        delete!(flat, matroid_element)
        term = find_flat_variable(chow_ring, flat)
        coloop_term_morphism = create_theta_i(coloop_term, chow_ring, matroid_element, term)
    end

    #TODO Change output into MChowRingDecomp
    first_term, second_term, coloop_term

end

"""
Creates theta_i as described on page 2 of
"A semi-small decomposition of the Chow ring of a matroid" by Braden et. al,
times a given factor in the image. E.g. x__1_3_5 * theta_i

These modified theta_i are the building blocks of the isomoprhisms of terms
inside the direct sum decomposition of a Chow ring with Chow Rings of smaller
matroids.
"""
function create_theta_i(domain::MChowRing, image::MChowRing, deleted_element::Int64, factor=1)
    # TODO: Replace "deleted element" with "deleted_set"!
    # The terms in the big sum are matroids which result from deleting many elements.

    image_gens = Vector{MPolyQuoElem{fmpq_mpoly}}(undef, length(gens(domain.chow_ring)))
    domain_gens = gens(domain.chow_ring)

    for i = 1:length(gens(domain.chow_ring))
        domain_gen = domain_gens[i]
        image_gens[i] = factor * find_gen_image(domain_gen, image, deleted_element)
    end

    hom(domain.chow_ring, image.chow_ring, image_gens)
end

# TODO: This works... but it's ugly string hacking.
# One should instead pass around a dict of set ints to variables instead.
"""
Helper function for create_theta_i. Given e.g. domain_gen = x_1_3_5_7 and
deleted element = 5, it creates the polynomial

x_1_3_6_8 + x_1_3_5_6_8 in the image, where a term is set to 0 if the
 corresponding variable does not exist in the image.
"""
function find_gen_image(domain_gen, image::MChowRing, deleted_element::Int64)
     elements_in_flat = split(string(domain_gen), "__")[2]

     first_canidate = "x_"
     second_canidate = "x_"

     previous_element = nothing
     second_canidate_adjusted = false
     for element_string in split(elements_in_flat, "_")
         element = parse(Int, element_string)
         if element >= deleted_element
             element += 1
         end
         if element > deleted_element &&
            (previous_element == nothing || previous_element < deleted_element)
              second_canidate = second_canidate * "_" * string(deleted_element)
              second_canidate_adjusted = true
         end

         first_canidate = first_canidate * "_" * string(element)
         second_canidate = second_canidate * "_" * string(element)
         previous_element = element
     end

     if !second_canidate_adjusted
         second_canidate = second_canidate * "_" * string(deleted_element)
     end
     first_term = image.chow_ring()
     second_term = image.chow_ring()
     for image_gen in gens(image.chow_ring)
         image_gen_string = string(image_gen)
         if image_gen_string == first_canidate
             first_term = image_gen
         end
         if image_gen_string == second_canidate
             second_term = image_gen
         end
     end

     first_term + second_term
end

# TODO: Find some way to remove this ugly string hacking
"""
Finds the generator in chow_ring with exactly the elements as in contents.
"""
function find_flat_variable(chow_ring::MChowRing, contents)
    for gen in gens(chow_ring.chow_ring)
        elements_in_flat = split(split(string(gen), "__")[2], "_")

        if length(elements_in_flat) != length(contents)
            continue
        end

        sets_equal = true
        for element_string in elements_in_flat
            element = parse(Int, element_string)
            if !Polymake.in(element, contents)
                sets_equal = false
                break
            end
        end

        if !sets_equal
            continue
        end

        # If the control flow reaches here, then the gen has the same elements
        # as contents.

        return gen
    end
end

function matroid_chow_ring(matroid::pm.BigObject)::MChowRing
    if pm.type_name(matroid) != "Matroid"
        throw(ArgumentError("BigObject is not a matroid."))
    end

    flats = matroid.LATTICE_OF_FLATS.FACES;
    #TODO: Is there a way to avoid this conversion?
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)

    # Technically non-empty proper flats. First element is either empty or the whole set, last element is the other of the two.
    proper_flats = flats[2:pm.size(flats, 1)-1, 1:pm.size(flats,2)]

    if length(proper_flats) == 0
        # TODO: Better representation of trivial ring?
        return MChowRing(nothing, [], matroid)
    end

    base_ring, indeterminates = generate_base_ring(proper_flats)

    type_i_polynomials = generate_type_i_ideal(base_ring, proper_flats, indeterminates, matroid)
    type_j_polynomials = generate_type_j_ideal(proper_flats, indeterminates)
    generators = vcat(type_i_polynomials, type_j_polynomials)

    chow_modulus = ideal(base_ring, generators)
    # Add std(ideal) to compute standard basis?
    chow_ring, projection = quo(base_ring, chow_modulus)
    projected_indeterminates = [projection(indeterminate) for indeterminate in indeterminates]

    MChowRing(chow_ring, projected_indeterminates, matroid)
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

    projected_element_vars = [projection(matroid_element_var) for matroid_element_var in matroid_element_vars]
    projected_flat_vars = [projection(flat_var) for flat_var in flat_vars]
    chow_ring, projected_element_vars, projected_flat_vars
end


function generate_base_ring(proper_flats)
    variable_names = create_flat_variables_names(proper_flats)
    # TODO: This doesn't seem to work well when it represents the trivial ring.
    #       Printing a trivial ring generated here causes errors.
    #       Perhaps a different way to create the trivial ring?
    PolynomialRing(QQ, variable_names);
end

function generate_augmented_base_ring(proper_flats, n_elements)
    matroid_element_variable_names = Array{String}(undef, n_elements)
    for i in 0:n_elements-1
        matroid_element_variable_names[i+1] = "y_" * string(i)
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
        row = extract_flat(flats, i)
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
                proper_flat = extract_flat(proper_flats, flat_index)
                print(proper_flat)
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
            flat = extract_flat(proper_flats, j)
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
            i_flat = extract_flat(proper_flats, i)
            j_flat = extract_flat(proper_flats, j)
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
        flat = extract_flat(proper_flats, i)
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

function extract_flat(flats, index)
    flat = Polymake.row(flats, index)

    # This is a hacky fix for the conversino to IncidenceMatrix adding 1 to every index.
    result = Array{Int64}(undef, length(flat))

    i = 1
    for element in flat
        result[i] = element - 1
        i += 1
    end

    Polymake.Set(result)
end

end
