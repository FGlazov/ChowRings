# TODO: Split this up into multiple files.

module MatroidChowRings

export direct_sum_decomp, matroid_chow_ring, augmented_matroid_chow_ring, apply_homorphism

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
    indeterminates#::Vector{MPolyQuoElem{fmpq_mpoly}}
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

first_term_morphism + coloop_term_morphism + sum_((projection_1, projection_2, term)) term * projection_1 * projection 2

Where g_f,h_f are projections from CH(M_F+1) and CH(M^F), respectively, to CH(M-i).
The input to g_f and h_f would then be values in the corresponding tuple in MChowRingDecomp.

Fields:
- first_term_morphism:  Sends elements of the first_term into the LHS.
- second_term_morphism: A vector of tuples of morphisms and a term x_(F+i), written
                        as (projection_1, projection_2, x_(F+i)). Each tuple represents
                        the psi described in Propositon 2.21 of Tom Braden et. al.
                        It is used in as isomorphism in the direct sum decompsotion, see
                        Propositon 3.5. It can be seen as the product of two projections, a theta_i
                        along with a constant term. I.e. as term * projection_1 * projection_2.
                        Here projection_2 already contains the theta_i.
- coloop_term_morphism: An optional morphism, that sends members of the
                        coloopterm into the LHS. Equal to X_(E-i) * theta_i
"""
struct MChowRingHomorphism
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
    homomorphism::MChowRingHomorphism #TODO: Rename this "isomorphism"?
end


"""
    function direct_sum_decomp(matroid::pm.BigObject, matroid_element::Int64)

Convinience function to compute the direct sum decomposition of the chow ring of a matroid without
first explicitely computing the chow ring itself. See documentation of the same function applied
to a chow ring of type MChowRing for more details.

The matroid may be any loopless matroid.
Matroid element is a number from 0 to matroid.NR_ELEMENTS-1, inclusive.
The matroid element is the one which is removed in the decomposition.

# Example
'''julia-repl
julia> using Oscar
julia> using MatroidChowRings
julia> f8 = Polymake.matroid.f8_matroid();
julia> decomp = direct_sum_decomp(f8, 3);
'''
"""
function direct_sum_decomp(matroid::pm.BigObject, matroid_element::Int64)
    direct_sum_decomp(matroid_chow_ring(matroid), matroid_element)
end

# TODO: Implement an augmented version of this.
"""
Computes the 


"""
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
                to_remove_1 = proper_flat_copy
                to_remove_2 = set_complement(ground_set, proper_flat)

                matroid_1 = pm.matroid.contraction(matroid, to_remove_1)
                matroid_2 = pm.matroid.deletion(matroid, to_remove_2)

                chow_1 = matroid_chow_ring(matroid_1)
                chow_2 = matroid_chow_ring(matroid_2)

                push!(second_term, (chow_1, chow_2))

                projection_1 = create_projection(chow_1, chow_ring, to_remove_1, matroid_element, true)
                projection_2 = create_projection(chow_2, chow_ring, to_remove_2, matroid_element, false)
                term = find_flat_variable(first_term, proper_flat_copy)
                push!(second_term_morphism, (projection_1, projection_2, term))
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

    chow_hom = MChowRingHomorphism(first_term_morphism, second_term_morphism, coloop_term_morphism)

    MChowRingDecomp(matroid_element, chow_ring, first_term, second_term, coloop_term, chow_hom)

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

"""
Creates a projection function from a chow ring to another chow ring, where the
corresponding domain matroid is a contraction or deletion of the codomain matroid.

This is a building block of the pushforward map psi described in Proposition 2.18
and used in Proposition 3.5 to describe the isomorphisms behind the Chow ring
decomposition.

For contraction: The removed elements correspond to a flat F+i. In this case,
X_(F') gets mapped to X_(F+i union F')

For deletion: The removed elements correspond to the complement of a flat F.
There still exists a corresponding i. In this case, X_(F') gets mapped to
X_(F') + X_(F'+i).
"""
function create_projection(domain::MChowRing, image::MChowRing, removed_elements, i::Int64, is_contraction::Bool)
    if domain.chow_ring == nothing # Edge case where ring is trivial.
        return nothing
    end

    image_gens = Vector{MPolyQuoElem{fmpq_mpoly}}(undef, length(gens(domain.chow_ring)))
    domain_gens = gens(domain.chow_ring)

    # Relates matroid elements of the domain with the image
    # This is the reindexing of contraction/deletion in reverse.
    int_mapping = Array{Int64}(undef, domain.matroid.N_ELEMENTS)
    next_free_element = 0
    while Polymake.in(next_free_element, removed_elements)
        next_free_element += 1
    end
    for j = 1:length(int_mapping)
        int_mapping[j] = next_free_element

        next_free_element += 1
        while Polymake.in(next_free_element, removed_elements)
            next_free_element += 1
        end
    end


    for j = 1:length(gens(domain.chow_ring))
        domain_gen = domain_gens[j]
        image_gens[j] = project_gen(domain_gen, image, removed_elements, i, is_contraction, int_mapping)
    end

    hom(domain.chow_ring, image.chow_ring, image_gens)
end

function project_gen(domain_gen::MPolyQuoElem{fmpq_mpoly}, image::MChowRing, removed_elements, i::Int64, is_contraction::Bool, int_mapping)
    elements_in_gen = split(split(string(domain_gen), "__")[2], "_")
    elements_in_gen = [parse(Int, e) for e in elements_in_gen]

    # +1 here because matroid elements start from 0, and arrays start from 1.
    contents_image = [int_mapping[e+1] for e in elements_in_gen]

    if is_contraction
        contents_image = vcat(contents_image, Vector(removed_elements))
        sort!(contents_image)
        return find_flat_variable(image, contents_image)
    else
        x1 = find_flat_variable(image, contents_image)
        append!(contents_image, i)
        sort!(contents_image)
        x2 = find_flat_variable(image, contents_image)
        # Due to the structure of S_i, x1 and x2 should always be found.
        return x1 + x2
    end

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

# TODO: Add types here. Test it. Add example.
"""
This applies the MChowRingHomorphism, sending an element of the RHS of a
direct sum decomposition into the LHS.

Note that the right hand side of the decomposition is not represented as
subalgebras of CH(M), but as Chow rings which are isomorphic to the subalgebras
of CH(M). If you wish to turn an element on the RHS to one which lives in the
domain of the LHS, you need to use this structure.

So using this function you can instead recover an element which truely does
belong to the original CH(M).

Note that the fields here are strongly correlated. The chow_ring_hom must be a
MChowRingHomorphism from some MChowRingDecomp.

The first_term field must be an element inside the first_term of MChowRingDecomp.
The second_term must be a vector of 2-tuples, containing elements of the
second_term of the corresponding second_term of a MChowRingDecomp, in the same
order as they are in MChowRingDecomp.

"""
function apply_homorphism(chow_ring_hom::MChowRingHomorphism, first_term, second_term, coloop_term)
    result = chow_ring_hom.first_term_morphism(first_term)

    i = 1
    for (contracted_term, deleted_term) in second_term
        projection_1, projection_2, term = chow_ring_hom.second_term_morphism[i]

        result +=  term * projection_1(contracted_term) * projection_2(deleted_term)
        i += 1
    end

    if coloop_term != nothing
        result += chow_ring_hom.coloop_term_morphism(coloop_term)
    end

    result
end

"""
    matroid_chow_ring(matroid::pm.BigObject)::MChowRing

Computes the Chow ring of a matroid. It follows the convention of the chow ring as defined
in 'A SEMI-SMALL DECOMPOSITION OF THE CHOW RING OF A MATROID' by Tom Braden, June Huh et. al.
The chow ring is represented over the rationals, and only nonempty proper flats enter as variables.

This function accepts any loopless matroid as input, where the matroid is a polymake matroid. 
It returns the Chow ring described as a quitotent ring, and it does not precompute more than
it needs to. In particular, if you do any calculations on the ring, then you will likely need
to wait a bit the first time as it will likely compute a Gröbner basis

# Example
'''julia-repl
julia> using Oscar
julia> using MatroidChowRings
julia> fano_matroid = Polymake.matroid.fano_matroid();
julia> chow_ring =  matroid_chow_ring(fano_matroid);
julia> x = chow_ring.indeterminates
14-element Vector{MPolyQuoElem{fmpq_mpoly}}:
 x__0
 x__1
 x__2
 x__3
 x__4
 x__5
 x__6
 x__0_3_4
 x__0_2_5
 x__0_1_6
 x__1_2_3
 x__1_4_5
 x__2_4_6
 x__3_5_6
 julia> x[1]
 x__0
 julia> x[8]
 x__0_3_4
 julia> x[11]
 x__1_2_3
 julia> x[1] * x[8]
 -x__3_5_6^2
 julia> x[1] * x[11]
 0

'''

"""
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

# TODO: Document and add example.
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
        set_name = split(string(extract_flat(flats, i)), "\n")[2]
        variable_descriptor = replace(set_name[2:length(set_name)-1], " " => "_")
        variable_names[i] = "x__" * variable_descriptor
    end

    variable_names
end


function generate_type_i_ideal(base_ring, proper_flats, indeterminates, matroid)
    n = matroid.N_ELEMENTS
    n_proper_flats = pm.size(proper_flats, 1)
    ideal_polynomials = Vector{fmpq_mpoly}()

    for i in 0:n-1
        for j in i+1:n-1
            ij_polynomial = base_ring() # Create zero.
            ij_set = Set{Int64}([i,j])
            for flat_index in 1:n_proper_flats
                proper_flat = extract_flat(proper_flats, flat_index)
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

    # This is a hacky fix for the conversion to IncidenceMatrix adding 1 to every index.
    result = Array{Int64}(undef, length(flat))

    i = 1
    for element in flat
        result[i] = element - 1
        i += 1
    end

    Polymake.Set(result)
end

end
