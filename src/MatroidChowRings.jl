module MatroidChowRings

using Oscar;
const pm = Polymake;


function matroid_chow_ring(matroid::pm.BigObject)
    if pm.type_name(matroid) != "Matroid"
        throw(ArgumentError("BigObject is not a matroid."))
    end

    flats = matroid.LATTICE_OF_FLATS.FACES;
    flats = pm.@pm common.convert_to{IncidenceMatrix}(flats)
    top_flat_index = matroid.LATTICE_OF_FLATS.TOP_NODE
    if top_flat_index == 0
        proper_flats = flats[2:pm.size(flats, 1), 1:pm.size(flats,2)]
    else
        proper_flats = flats[1:pm.size(flats, 1)-1, 1:pm.size(flats,2)]
    end

    R, indeterminates = generate_base_ring(proper_flats)
    generate_type_i_ideal(proper_flats, indeterminates, matroid)


    print("Hello chow rings!")
end

function generate_base_ring(proper_flats)
    #variable_names = ["x__" * replace(split(string(proper_flat), "\n")[2], " " => "_")
    #                    for proper_flat in proper_flats]
    n_proper_flats = pm.size(proper_flats, 1)
    variable_names = Array{String}(undef, n_proper_flats)

    for i in 1:n_proper_flats
        row = pm.row(proper_flats, i)
        # TODO: Is there some better way to extract the polymake string without the type?
        set_name = split(string(pm.row(proper_flats, i)), "\n")[2]
        variable_descriptor = replace(set_name[2:length(set_name)-1], " " => "_")
        variable_names[i] = "x__" * variable_descriptor
    end

    Singular.PolynomialRing( Singular.QQ, variable_names);
end

function generate_type_i_ideal(proper_flats, indeterminates, matroid)
    n = matroid.N_ELEMENTS
    n_proper_flats = pm.size(proper_flats, 1)
    ideal_polynomials = Vector{Singular.spoly{Singular.n_Q}}()
    for i in 1:n
        for j in i+1:n
            ij_polynomial = indeterminates[1] - indeterminates[1] # TODO: Better way to create zero polynomial?
            ij_set = Set{Int64}([i,j])
            for flat_index in 1:n_proper_flats
                proper_flat = pm.row(proper_flats, flat_index)
                if pm.in(i, proper_flat)
                    ij_polynomial += indeterminates[flat_index]
                end
                if pm.in(i, proper_flat)
                    ij_polynomial -= indeterminates[flat_index]
                end
            end

            append!(ideal_polynomials, ij_polynomial)
        end
    end



end


end
