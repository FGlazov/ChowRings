# TODO: Split this up into multiple files.

module MatroidChowRings

export direct_sum_decomp, matroid_chow_ring, apply_homorphism
export augmented_direct_sum_decomp, augmented_matroid_chow_ring


using Oscar;
const pm = Polymake;

"""
    struct MChowRing

Represents a matroid Chow ring. The result of a matroid_chow_ring call.

Fields
- chow_ring:      The Chow ring of the matroid, represented as an affine algebra.
                  It may be == nothing, in that case chow ring is QQ.
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
    indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}
    matroid::pm.BigObject
    chow_ring::MPolyQuo{MPolyElem_dec{fmpq, fmpq_mpoly}}

    function MChowRing(indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}, matroid::pm.BigObject)
        if pm.type_name(matroid) != "Matroid"
            error("BigObject is not a matroid")
        end

        new(indeterminates, matroid)
    end


    function MChowRing(indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}, matroid::pm.BigObject, chow_ring)
        if pm.type_name(matroid) != "Matroid"
            error("BigObject is not a matroid")
        end
        new(indeterminates, matroid, chow_ring)
    end
end

"""
    struct MAugChowRing

Represents am augmented matroid Chow ring. The result of a augmented_matroid_chow_ring call.

Fields
- chow_ring:              The augmented Chow ring of the matroid, represented as an affine algebra.
                          It may be == nothing, in that case chow ring is QQ.
- flat_indeterminates:    The flat indeterminates of the augmented Chow ring. These are of the form
                          x_F, where F is a proper flat. The name of F corresponds to elements inside
                           that flat. E.g. x__1_5_7 would be the flat consisting of the
                          elements 1, 5, and 7. Another way to see all the flats would be to call
                          matroid.LATTICE_OF_FLATS.FACES.
                          This vector may be empty, in that case the chow ring is QQ.
- element_indeterminates: The indeterminates of the augmented Chow ring. These are of the form y_i,
                          where i is an element of the groundset of the matroid. The number of these
                          corresponds to matroid.N_ELEMENTS, ranging from 0 to N_ELEMENTS-1.
- matroid:                The matroid whose augmented Chow ring has been computed. It is
                          represented as a polymake Matroid.
"""
struct MAugChowRing
    flat_indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}
    element_indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}
    matroid::pm.BigObject
    chow_ring::MPolyQuo{MPolyElem_dec{fmpq, fmpq_mpoly}}

    function MAugChowRing(flat_indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}, element_indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}, matroid::pm.BigObject)
        if pm.type_name(matroid) != "Matroid"
            error("BigObject is not a matroid")
        end

        new(flat_indeterminates, element_indeterminates, matroid)
    end

    function MAugChowRing(flat_indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}, element_indeterminates::Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}, matroid::pm.BigObject, chow_ring)
        if pm.type_name(matroid) != "Matroid"
            error("BigObject is not a matroid")
        end
        new(flat_indeterminates, element_indeterminates, matroid, chow_ring)
    end
end

"""
    struct MChowRingHomorphism

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
    struct MAugChowRingHomorphism

Represents a homorphism which sends the RHS of an augmentedChow ring direct sum
decomposition to the LHS of the decompositon.

Note that the right hand side of the decomposition is not represented as
subalgebras of CH_aug(M), but as (augmented) Chow rings which are isomorphic to the
subalgebras of CH_aug(M). If you wish to turn an element on the RHS to one
which lives in the domain of the LHS, you need to use this structure.

This structure represents the morphism:

first_term_morphism + coloop_term_morphism + sum_((projection_1, projection_2, term)) term * projection_1 * projection 2

Where g_f,h_f are projections from CH(M_F+1) and CH_aug(M^F), respectively, to CH_aug(M-i).
The input to g_f and h_f would then be values in the corresponding tuple in MAugChowRingDecomp.

Fields:
- first_term_morphism:  Sends elements of the first_term into the LHS.
- second_term_morphism: A vector of tuples of morphisms and a term x_(F+i), written
                        as (projection_1, projection_2, x_(F+i)). Each tuple represents
                        the psi described in Propositon 2.24 of Tom Braden et. al.
                        It is used in as isomorphism in the direct sum decompsotion, see
                        Propositon 3.5. It can be seen as the product of two projections, a theta_i
                        along with a constant term. I.e. as term * projection_1 * projection_2.
                        Here projection_2 already contains the theta_i.
- coloop_term_morphism: An optional morphism, that sends members of the
                        coloopterm into the LHS. Equal to X_(E-i) * theta_i
"""
struct MAugChowRingHomorphism
    first_term_morphism
    second_term_morphism
    coloop_term_morphism
end

"""
    struct MChowRingDecomp
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
                   i is the deleted element, and F is a flat. This is the part
                   of the decomposition indexed by the S_i in the paper.
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

"""
    struct MAugChowRingDecomp
Represents a direct sum decomposition of the augmented Chow ring of a matroid. See
"A semi-small decomposition of the Chow ring of a matroid by Braden et. al" for
details. The result of a direct_sum_decomp call.

Note that the right hand side of the decomposition is not represented as
subalgebras of CH_aug(M), but as Chow rings which are isomorphic to the subalgebras
of CH_aug(M). If you wish to turn an element on the RHS to one which lives in the
domain of the LHS, you need to use the homomorphism field.

Fields:
- deleted_element: The element which was deleted on the RHS of the decomposition.
- chow_ring_LHS:   The augmented Chow ring which was decomposed.
- first_term:      Always present part of the decomposition, equal to the augmented Chow
                   ring of M-i.
- second_term:     A vector of tuples of chow rings. This has varying length,
                   and each tuple is of the form (CH(M_(F+i)), CH_aug(M^F), where
                   i is the deleted element, and F is a flat. This is the part
                   of the decomposition indexed by the S_i in the paper.
- coloopterm:      A term that appears in the decomposition if deleted_element
                   is a coloop. Equal to the chow ring of M-i if present.
- homomorphism:    Represents the homomorphism which takes a combination of
                   elements from the first/second_term/coloop_term and maps it
                   to the LHS of the decomposition.
"""
struct MAugChowRingDecomp #TODO Add types here.
    deleted_element::Int64
    chow_ring_LHS::MAugChowRing
    first_term::MAugChowRing
    second_term
    coloopterm
    homomorphism::MAugChowRingHomorphism
end


"""
    function direct_sum_decomp(matroid::pm.BigObject, matroid_element::Int64)

Convinience function to compute the direct sum decomposition of the Chow ring of a matroid without
first explicitely computing the chow ring itself. See documentation of the same function applied
to a chow ring of type MChowRing for more details.

The matroid may be any loopless matroid.
Matroid element is a number from 0 to matroid.NR_ELEMENTS-1, inclusive.
The matroid element is the one which is removed in the decomposition.

# Example
```julia-repl
julia> using Oscar
julia> using MatroidChowRings
julia> f8 = Polymake.matroid.f8_matroid();
julia> decomp = direct_sum_decomp(f8, 3);
julia> length(gens(decomp.chow_ring_LHS.chow_ring))
56
julia> length(gens(decomp.first_term.chow_ring))
45
julia> length(decomp.second_term)
10
julia> decomp.coloopterm == nothing
true
```

Note that the length of the gens corresponds to the number of nonempty proper flats.
The 56 and 45 above say that the matroid corresponding to f8 delete 3 has 13 less flats.

Each element of decomp.second_term corresponds to a flat F, such that F+3 is also a flat.
Thus the length 10 says that there are 10 such choices for F.

f8 has no coloops, so the coloopterm must be empty.
"""
function direct_sum_decomp(matroid::pm.BigObject, matroid_element::Int64)
    direct_sum_decomp(matroid_chow_ring(matroid), matroid_element)
end

# TODO: Implement an augmented version of this.
"""
    function direct_sum_decomp(chow_ring::MChowRing, matroid_element::Int64)

Computes the (semi-small) direct sum decomposition of the Chow ring of a matroid, as first described
in 'A SEMI-SMALL DECOMPOSITION OF THE CHOW RING OF A MATROID' by Tom Braden, June Huh et. al.

Due to technical reasons the decomposition is not represented as
subalgebras of CH(M), but as Chow rings which are isomorphic to the subalgebras
of CH(M). Meaning that instead of a direct sum decomposition, this function instead
returns all the (isomorphic) components of the decomposition.

If you wish to instead recover the associated subalgebras, the output of this function also
contains the relevant pushforward map, which takes an element from the product of all the
isomoprhic components of the decomposition, and maps it to an element of the matroid chow ring
in the image. By then mapping a set of generators in the isomorphic component, you can get the
corresponding generators of the subalgebra. See the function apply_homorphism for details.

Note that because polymake reindexes matroid elements to always run from 0 to n-1, an indeterminate
of the form e.g x_1_3_5 need not correspond to x_1_3_5 in the original matroid. In particular, if
3 was removed in the construction of the smaller matroid. The function apply_homorphism accounts for
the reindexing.

chow_ring         the result of a direct_sum_decomp call
matroid_element   the matroid element to be removed in the decomposition.

# Example

```julia-repl
julia> using Oscar
julia> using MatroidChowRings
julia> fano_matroid = Polymake.matroid.fano_matroid();
julia> chow_ring =  matroid_chow_ring(fano_matroid);
julia> x = chow_ring.indeterminates;
julia> decomp = direct_sum_decomp(chow_ring, 0);
julia> length(decomp.first_term.indeterminates)
13
julia> length(decomp.second_term)
0
julia> decomp.coloopterm == nothing
true
julia> morphism = decomp.homomorphism.first_term_morphism;
julia> isbijective(morphism)
true
julia> morphism.image
13-element Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}:
 x__6 - x__1_2_3 - x__1_4_5 + x__2_4_6 + x__3_5_6
 x__6 - x__0_2_5 + x__0_1_6 - x__1_2_3 + x__3_5_6
 x__6 - x__0_3_4 + x__0_1_6 - x__1_2_3 + x__2_4_6
 x__6 - x__0_3_4 + x__0_1_6 - x__1_4_5 + x__3_5_6
 x__6 - x__0_2_5 + x__0_1_6 - x__1_4_5 + x__2_4_6
 x__6
 x__1_2_3
 x__1_4_5
 x__0_1_6
 x__0_2_5
 x__2_4_6
 x__0_3_4
 x__3_5_6
 julia> for img_gen in morphism.image
           if img_gen in x
               print(true,"\n")
           else
               print(false,"\n")
           end
       end
true
true
true
true
true
true
true
true
true
true
true
true
true
```

The above code says that the direct sum decomposition of the fano matroid consists of only
one component which is isomorphic to matroid you get by deleting any element from the fano
matroid. This means that the pushforward map is very simple - it is composed of a single
algebraic homomorphism. One can directly verify that this is indeed a decomposition by
checking that this algerbraic homoprhism is isomorphic.

The for loop proves that generators in the RHS of the decomposition get mapped to generators
in the Chow ring of the Fano matroid. This paired with the fact that the morphism is isomorphic,
shows that the natural set of generators for the Fano matroid contain redudant elements, in particular
there is a generating set of length 13.
"""
function direct_sum_decomp(chow_ring::MChowRing, matroid_element::Int64)
    if matroid_element < 0 || matroid_element >= chow_ring.matroid.N_ELEMENTS
        print("Matroid element out of bounds, can not decompose.\n")
        return nothing
    end
    matroid = chow_ring.matroid

    first_term = matroid_chow_ring(pm.matroid.deletion(matroid, matroid_element))
    first_term_morphism = create_theta_i(first_term, chow_ring, matroid_element)

    ground_set = Polymake.Set(range(0, matroid.N_ELEMENTS - 1, step=1))
    flats = matroid.LATTICE_OF_FLATS.FACES
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)
    proper_flats = flats[2:pm.size(flats, 1)-1, 1:pm.size(flats,2)]
    n_proper_flats = pm.size(proper_flats, 1)

    flat_canidates = []
    for flat_index in true_indices_in_col(proper_flats, matroid_element + 1)
        proper_flat = extract_flat(proper_flats, flat_index)
        push!(flat_canidates, proper_flat)
    end
    flat_canidates = Set(flat_canidates)
    second_term = []
    second_term_morphism = []

    for flat_index in false_indices_in_col(proper_flats, matroid_element + 1)
        proper_flat = extract_flat(proper_flats, flat_index)
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
            term = find_flat_variable(chow_ring, proper_flat_copy)
            push!(second_term_morphism, (projection_1, projection_2, term))
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

function augmented_direct_sum_decomp(matroid::pm.BigObject, matroid_element::Int64)
    augmented_direct_sum_decomp(augmented_matroid_chow_ring(matroid), matroid_element)
end

function augmented_direct_sum_decomp(chow_ring::MAugChowRing, matroid_element::Int64)
    if matroid_element < 0 || matroid_element >= chow_ring.matroid.N_ELEMENTS
        print("Matroid element out of bounds, can not decompose.\n")
        return nothing
    end

    matroid = chow_ring.matroid
    first_term = augmented_matroid_chow_ring(pm.matroid.deletion(matroid, matroid_element))
    first_term_morphism = create_aug_theta_i(first_term, chow_ring, matroid_element)

    ground_set = Polymake.Set(range(0, matroid.N_ELEMENTS - 1, step=1))
    flats = matroid.LATTICE_OF_FLATS.FACES
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)

    top_node = matroid.LATTICE_OF_FLATS.TOP_NODE
    if top_node == 0
        proper_flats = flats[2:pm.size(flats, 1), 1:pm.size(flats,2)]
    else
        proper_flats = flats[1:pm.size(flats, 1)-1, 1:pm.size(flats,2)]
    end

    flat_canidates = []
    for flat_index in true_indices_in_col(proper_flats, matroid_element + 1)
        proper_flat = extract_flat(proper_flats, flat_index)
        push!(flat_canidates, proper_flat)
    end

    second_term = []
    second_term_morphism = []

    for flat_index in false_indices_in_col(proper_flats, matroid_element + 1)
        proper_flat = extract_flat(proper_flats, flat_index)
        proper_flat_copy = copy(proper_flat)
        push!(proper_flat_copy, matroid_element)
        if proper_flat_copy in flat_canidates
            to_remove_1 = proper_flat_copy
            to_remove_2 = set_complement(ground_set, proper_flat)

            matroid_1 = pm.matroid.contraction(matroid, to_remove_1)
            matroid_2 = pm.matroid.deletion(matroid, to_remove_2)

            chow_1 = matroid_chow_ring(matroid_1)
            chow_2 = augmented_matroid_chow_ring(matroid_2)

            push!(second_term, (chow_1, chow_2))

            projection_1 = create_projection(chow_1, chow_ring, to_remove_1, matroid_element, true)
            projection_2 = create_projection(chow_2, chow_ring, to_remove_2, matroid_element, false)
            term = find_flat_variable(chow_ring, proper_flat_copy)
            push!(second_term_morphism, (projection_1, projection_2, term))
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
        coloop_term_morphism = create_aug_theta_i(coloop_term, chow_ring, matroid_element, term)
    end

    chow_hom = MAugChowRingHomorphism(first_term_morphism, second_term_morphism, coloop_term_morphism)

    MAugChowRingDecomp(matroid_element, chow_ring, first_term, second_term, coloop_term, chow_hom)
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
    image_gens = Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}(undef, length(gens(domain.chow_ring)))
    domain_gens = gens(domain.chow_ring)

    for i = 1:length(gens(domain.chow_ring))
        domain_gen = domain_gens[i]
        image_gens[i] = factor * find_gen_image(domain_gen, image, deleted_element)
    end

    hom(domain.chow_ring, image.chow_ring, image_gens)
end

function create_aug_theta_i(domain::MAugChowRing, image::MAugChowRing, deleted_element::Int64, factor=1)
    image_gens = Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}(undef, length(gens(domain.chow_ring)))
    domain_ys = domain.element_indeterminates
    domain_xs = domain.flat_indeterminates
    image_ys = image.element_indeterminates

    for i = 1:length(domain_ys)
        image_index = i
        if i > deleted_element
            image_index += 1
        end
        image_gens[i] = factor * image_ys[image_index]
    end

    offset = length(domain_ys)
    for i = 1:length(domain_xs)
        domain_gen = domain_xs[i]
        image_gens[i + offset] = factor * find_gen_image(domain_gen, image, deleted_element)
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
function find_gen_image(domain_gen, image, deleted_element::Int64)
     elements_in_flat = split(string(domain_gen), "__")[2]

     first_canidate = "x_"
     second_canidate = "x_"

     previous_element = nothing
     second_canidate_adjusted = false

     for element_string in split(elements_in_flat, "_")
         if element_string == ""
             break
         end

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
    if !isdefined(domain, :chow_ring) # Edge case where ring is trivial.
        return nothing
    end

    image_gens = Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}(undef, length(gens(domain.chow_ring)))
    domain_gens = gens(domain.chow_ring)
    int_mapping = create_int_mapping(domain.matroid.N_ELEMENTS, removed_elements)

    for j = 1:length(gens(domain.chow_ring))
        domain_gen = domain_gens[j]
        image_gens[j] = project_gen(domain_gen, image, removed_elements, i, is_contraction, int_mapping)
    end

    hom(domain.chow_ring, image.chow_ring, image_gens)
end

function create_projection(domain::MAugChowRing, image::MAugChowRing, removed_elements, i::Int64, is_contraction::Bool)
    if !isdefined(domain, :chow_ring) # Edge case where ring is trivial.
        return nothing
    end

    image_gens = Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}(undef, length(gens(domain.chow_ring)))
    domain_gens = gens(domain.chow_ring)
    int_mapping = create_int_mapping(domain.matroid.N_ELEMENTS, removed_elements)

    for j = 1:length(domain.element_indeterminates)
        domain_gen = domain.element_indeterminates[j]
        image_gens[j] = image.element_indeterminates[int_mapping[j] + 1]
    end

    offset = length(domain.element_indeterminates)
    for j = 1:length(domain.flat_indeterminates)
        domain_gen = domain.flat_indeterminates[j]
        image_gens[j + offset] = project_gen(domain_gen, image, removed_elements, i, is_contraction, int_mapping)
    end

    hom(domain.chow_ring, image.chow_ring, image_gens)
end

function create_projection(domain::MChowRing, image::MAugChowRing, removed_elements, i::Int64, is_contraction::Bool)
    if !isdefined(domain, :chow_ring) # Edge case where ring is trivial.
        return nothing
    end

    image_gens = Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}(undef, length(gens(domain.chow_ring)))
    domain_gens = gens(domain.chow_ring)
    int_mapping = create_int_mapping(domain.matroid.N_ELEMENTS, removed_elements)

    for j = 1:length(gens(domain.chow_ring))
        domain_gen = domain_gens[j]
        image_gens[j] = project_gen(domain_gen, image, removed_elements, i, is_contraction, int_mapping)
    end

    hom(domain.chow_ring, image.chow_ring, image_gens)
end

function create_int_mapping(n_elements::Int64, removed_elements)
    # Relates matroid elements of the domain with the image
    # This is the reindexing of contraction/deletion in reverse.
    int_mapping = Array{Int64}(undef, n_elements)
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

    int_mapping
end

function project_gen(domain_gen::MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}, image, removed_elements, i::Int64, is_contraction::Bool, int_mapping)
    elements_in_gen = split(split(string(domain_gen), "__")[2], "_")
    if elements_in_gen[1] == ""
        elements_in_gen = []
    end

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
function find_flat_variable(chow_ring, contents)
    for gen in gens(chow_ring.chow_ring)
        if split(string(gen), "__")[1] != "x"
            continue
        end

        elements_in_flat = split(split(string(gen), "__")[2], "_")
        if elements_in_flat[1] == ""
            elements_in_flat = []
        end

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

        if sets_equal
            return gen
        end
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

So using this function you can instead recover an element which truly does
belong to the original CH(M).

chow_ring_hom must be the homorphism which comes as one of the return values
of direct_sum_decomp.

The first_term field must be an element inside the first_term of MChowRingDecomp.
The second_term must be a vector of 2-tuples, containing elements of the
second_term of the corresponding second_term of a MChowRingDecomp, in the same
order as they are in MChowRingDecomp.

# Example

```julia-repl
julia> using Oscar
julia> using MatroidChowRings
julia> u45 = Polymake.matroid.uniform_matroid(4,5);
julia> decomp = direct_sum_decomp(u45, 2);
julia> length(decomp.second_term)
10
julia> decomp.coloopterm == nothing
true
julia> p1 = decomp.first_term.indeterminates[2] + decomp.first_term.indeterminates[1]
x__0 + x__1
julia> p2s = [];
julia> for (t1, t2) in decomp.second_term
           x1 = 1
           x2 = 1
           if length(t1.indeterminates) > 0
               x1 = t1.indeterminates[1]
           end
           if length(t2.indeterminates) > 0
               x2 = t2.indeterminates[2]
           end

           push!(p2s, (x1, x2))
       end
julia> p2s
10-element Vector{Any}:
(x__0, 1)
(x__0, 1)
(x__0, 1)
(x__0, 1)
(1, x__1)
(1, x__1)
(1, x__1)
(1, x__1)
(1, x__1)
(1, x__1)
julia> p3 = nothing;
julia> apply_homorphism(decomp.homomorphism, p1, p2s, p3)
2*x__4 - x__0_1*x__0_1_4 - 2*x__0_1 + x__0_2*x__0_2_4 - x__0_3*x__0_3_4 - x__0_3 - x__0_4*x__0_3_4 + x__0_4 + x__1_2*x__1_2_4 - x__1_3*x__1_3_4 - x__1_3 - x__1_4*x__1_3_4 + x__1_4 + x__2_3*x__2_3_4 + x__2_4*x__2_3_4 + 2*x__2_4 - x__3_4*x__2_3_4 + 2*x__3_4 - x__0_1_2^2 - 2*x__0_1_2 - 2*x__0_1_3 - x__0_2_3^2 - x__0_2_3 - x__0_2_4^2 + x__0_2_4 + x__0_3_4 - x__1_2_3^2 - x__1_2_3 - x__1_2_4^2 + x__1_2_4 + x__1_3_4 - x__2_3_4^2 + 2*x__2_3_4
```

"""
function apply_homorphism(chow_ring_hom::MChowRingHomorphism, first_term, second_term, coloop_term)
    result = chow_ring_hom.first_term_morphism(first_term)

    i = 1
    for (contracted_term, deleted_term) in second_term
        projection_1, projection_2, term = chow_ring_hom.second_term_morphism[i]

        mapped_contracted_term = contracted_term
        mapped_deleted_term = deleted_term
        if projection_1 != nothing
            mapped_contracted_term = projection_1(contracted_term)
        end
        if projection_2 != nothing
            mapped_deleted_term = projection_2(deleted_term)
        end

        result +=  term * mapped_contracted_term * mapped_deleted_term
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
```julia-repl
julia> using Oscar
julia> using MatroidChowRings
julia> fano_matroid = Polymake.matroid.fano_matroid();
julia> chow_ring =  matroid_chow_ring(fano_matroid);
julia> x = chow_ring.indeterminates
14-element Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}:
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
```

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
        empty_vector = Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}()
        return MChowRing(empty_vector, matroid)
    end

    base_ring, indeterminates = generate_base_ring(proper_flats)

    type_i_polynomials = generate_type_i_ideal(base_ring, proper_flats, indeterminates, matroid)
    type_j_polynomials = generate_type_j_ideal(proper_flats, indeterminates)
    generators = vcat(type_i_polynomials, type_j_polynomials)

    chow_modulus = ideal(base_ring, generators)
    # Add std(ideal) to compute standard basis?
    chow_ring, projection = quo(base_ring, chow_modulus)
    projected_indeterminates = [projection(indeterminate) for indeterminate in indeterminates]

    MChowRing(projected_indeterminates, matroid, chow_ring)
end

"""
    augmented_matroid_chow_ring(matroid::pm.BigObject)::MAugChowRing

Computes the augmented Chow ring of a matroid. It follows the convention of the chow ring as defined
in 'A SEMI-SMALL DECOMPOSITION OF THE CHOW RING OF A MATROID' by Tom Braden, June Huh et. al.
The chowVector ring is represented over the rationals, and proper flats and matroid groundset elements
enter as variables.

This function accepts any loopless matroid as input, where the matroid is a polymake matroid.
It returns the Chow ring described as a quitotent ring, and it does not precompute more than
it needs to. In particular, if you do any calculations on the ring, then you will likely need
to wait a bit the first time as it will likely compute a Gröbner basis.

# Example
```julia-repl
julia> using Oscar
julia> using MatroidChowRings
julia> fano_matroid = Polymake.matroid.fano_matroid();
julia> chow_ring = augmented_matroid_chow_ring(fano_matroid);
julia> x = chow_ring.flat_indeterminates
15-element Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}:
 x__
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

julia> y = chow_ring.element_indeterminates
7-element Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}:
 y_0
 y_1
 y_2
 y_3
 y_4
 y_5
 y_6
julia> y[1] * y[1]
0
julia> y[1] * y[2]
-x__6^2 - 4*x__6*x__3_5_6 - x__2_4_6^2 - x__3_5_6^2
julia> x[9]
x__0_3_4
julia> y[1] * x[9]
-x__0*x__0_1_6 - x__0_3_4^2
julia> y[2] * x[9]
0
```

"""
function augmented_matroid_chow_ring(matroid::pm.BigObject)::MAugChowRing
    if pm.type_name(matroid) != "Matroid"
        throw(ArgumentError("BigObject is not a matroid."))
    end

    if matroid.N_ELEMENTS == 0
        empty_vector = Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}()
        return MAugChowRing(empty_vector, empty_vector, matroid)
    end

    n_elements = matroid.N_ELEMENTS
    flats = matroid.LATTICE_OF_FLATS.FACES
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)

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
    chow_ring, projection = quo(base_ring, chow_modulus)

    projected_element_vars = [projection(matroid_element_var) for matroid_element_var in matroid_element_vars]
    projected_flat_vars = [projection(flat_var) for flat_var in flat_vars]

    MAugChowRing(projected_flat_vars, projected_element_vars, matroid, chow_ring)
end


function generate_base_ring(proper_flats)
    variable_names = create_flat_variables_names(proper_flats)
    GradedPolynomialRing(QQ, variable_names);
end

function generate_augmented_base_ring(proper_flats, n_elements)
    matroid_element_variable_names = Array{String}(undef, n_elements)
    for i in 0:n_elements-1
        matroid_element_variable_names[i+1] = "y_" * string(i)
    end

    flat_variable_names = create_flat_variables_names(proper_flats)
    variable_names = vcat(matroid_element_variable_names, flat_variable_names)

    ring, variables = GradedPolynomialRing( QQ, variable_names);

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
    ideal_polynomials = Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}()

    for i in 0:n-1
        for j in i+1:n-1
            ij_polynomial = base_ring() # Create zero.
            ij_set = Set{Int64}([i,j])

            for i_index in true_indices_in_col(proper_flats, i + 1)
                ij_polynomial += indeterminates[i_index]
            end

            for j_index in true_indices_in_col(proper_flats, j + 1)
                ij_polynomial -= indeterminates[j_index]
            end
            push!(ideal_polynomials, ij_polynomial)
        end
    end

    ideal_polynomials
end

function generate_augmented_type_i_ideal(proper_flats, matroid_element_vars, flat_vars)
    n_proper_flats = pm.size(proper_flats, 1)
    ideal_polynomials = Array{MPolyElem_dec{fmpq, fmpq_mpoly}}(undef, length(matroid_element_vars))

    ground_set = Polymake.Set(range(1, length(matroid_element_vars), step=1))


    for i in 1:length(matroid_element_vars)
        ideal_polynomials[i] = matroid_element_vars[i]

        for flat_index_no_i in false_indices_in_col(proper_flats, i)
            ideal_polynomials[i] -= flat_vars[flat_index_no_i]
        end
    end

    ideal_polynomials
end

function generate_type_j_ideal(proper_flats, indeterminates)
    n_proper_flats = pm.size(proper_flats, 1)
    ideal_polynomials = Vector{}()

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
    xy_polynomials = Vector{}()

    for i in 1:length(matroid_element_vars)
        for j in false_indices_in_col(proper_flats, i)
            push!(xy_polynomials, matroid_element_vars[i] * flat_vars[j])
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

function true_indices_in_col(matrix, index)
     Polymake.SparseArrays.nonzeroinds(matrix[:,index])
end

function false_indices_in_col(matrix, index)
     to_complement = true_indices_in_col(matrix, index)
     ground_set = Polymake.Set(range(1, Polymake.nrows(matrix), step=1))
     set_complement(ground_set, to_complement)
end

end
