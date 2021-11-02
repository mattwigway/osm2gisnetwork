# create turn features
# schema:
# https://desktop.arcgis.com/en/arcmap/latest/extensions/network-analyst/turns-in-the-network-dataset.htm

using OSMPBF
using Logging
using ProgressMeter
using Graphs
using MetaGraphs
using DataStructures
using Printf

# a turn is a left turn if between these headings
# note that these do overlap, because there is some ambiguity in whether some turns
# are straight/right/left etc, and may depend on perceptions.
const LEFT_TURN_RANGES = [(-170, -15)]
const STRAIGHT_RANGES = [(-30, 30)]
const RIGHT_TURN_RANGES = [(15, 170)]
const U_TURN_RANGES = [(-Inf32, -95), (95, Inf32)]

struct TurnRestriction
    segments::AbstractVector{EdgeRef}
    back::AbstractVector{Bool}
    geom::Vector{LatLon} 
    turn_angle::Float32
    osm_id::Int64
end

# Get the restriction type (no_left_turn etc) for a restriction
function get_rtype(r)
    if haskey(r.tags, "restriction")
        return r.tags["restriction"]
    elseif haskey(r.tags, "restriction:motorcar")
        return r.tags["restriction:motorcar"]
    else
        error("unreachable code")
    end
end

function write_turn_features(infile::String, outfile::String, way_segment_idx::Dict{Int64, Vector{EdgeRef}},
        node_locations::Dict{Int64, LatLon}, edges_for_node::Dict{Int64, Vector{EdgeRef}})

    n_wildcard = 0
    turn_restrictions = Vector{TurnRestriction}()
    @info "Parsing turn restrictions"
    rprog = ProgressUnknown()
    scan_pbf(infile, relations=r -> begin
        # first, check if this is a turn restriction
        if haskey(r.tags, "type") && r.tags["type"] == "restriction" && any(haskey.([r.tags], ["restriction", "restriction:motorcar"]))
            ProgressMeter.next!(rprog)

            rtype = get_rtype(r)

            if rtype == "no_exit" || rtype == "no_entry"
                n_wildcard += 1
            else
                # parse the restriction
                from = nothing
                to = nothing
                via = Vector{OSMPBF.RelationMember}()

                # sort members into roles
                for member in r.members
                    if member.role == "from"
                        if !isnothing(from)
                            @warn "Relation $(r.id) has multiple from members"
                        end
                        from = member
                    elseif member.role == "to"
                        if !isnothing(to)
                            @warn "Relation $(r.id) has multiple to members"
                        end
                        to = member
                    elseif member.role == "via"
                        push!(via, member)
                    else
                        @warn "Ignoring member with role $(member.role) in relation $(r.id)"
                    end
                end

                if isnothing(from)
                    @warn "Relation $(r.id) is missing from way"
                    return
                end

                if isnothing(to)
                    @warn "Relation $(r.id) is missing to way"
                    return
                end

                parsed = nothing
                if length(via) == 1 && via[1].type == OSMPBF.node
                    parsed = process_simple_restriction(r, from, to, via[1], way_segment_idx, node_locations)
                elseif length(via) ≥ 1 && all(map(v -> v.type == OSMPBF.way, via))
                    parsed = process_complex_restriction(r, from, to, way_segment_idx, node_locations)
                elseif length(via) == 0
                    @warn "restriction $(r.id) has no via members, skipping"
                else
                    @warn "via members of restriction $(r.id) are invalid (multiple nodes, mixed nodes/ways, relation members), skipping"
                end

                if !isnothing(parsed)
                    if startswith(rtype, "no_")
                        push!(turn_restrictions, parsed)
                    elseif startswith(rtype, "only_")
                        restrictions = convert_restriction_to_only_turn(parsed, node_locations, edges_for_node)
                        append!(turn_restrictions, restrictions)
                    else
                        @error "Skipping turn restriction of type $rtype"
                    end
                end

            end
        end
    end)

    @info "created $(length(turn_restrictions)) turn restrictions"

    if n_wildcard > 0
        @warn "Ignored $n_wildcard no_entry or no_exit restrictions"
    end

    # find the maximum number of edges in turn restrictions
    max_edges = max(map(r -> length(r.segments), turn_restrictions)...)

    @info "all turn restrictions have $max_edges edges or less"

    # https://discourse.julialang.org/t/how-to-create-a-new-shapefile-containing-a-few-points/43454/3
    ArchGDAL.create(outfile, driver = ArchGDAL.getdriver("ESRI Shapefile")) do ds
        # EPSG 4326 - WGS 84 - coordinate reference system used by OpenStreetMap
        ArchGDAL.createlayer(geom=ArchGDAL.wkbLineString, spatialref=ArchGDAL.importEPSG(4326)) do layer
            ArchGDAL.addfielddefn!(layer, "ObjectID", ArchGDAL.OFTInteger)
            ArchGDAL.addfielddefn!(layer, "Bearing", ArchGDAL.OFTReal)
            ArchGDAL.addfielddefn!(layer, "restricted", ArchGDAL.OFTInteger)
            ArchGDAL.addfielddefn!(layer, "OSMID", ArchGDAL.OFTInteger64)
            #ArchGDAL.addfielddefn!(layer, "Edge1End", ArchGDAL.OFTInteger)
            # for i in 1:max_edges
            #     ArchGDAL.addfielddefn!(layer, "Edge$(i)FID", ArchGDAL.OFTInteger)
            #     ArchGDAL.addfielddefn!(layer, "Edge$(i)FCID", ArchGDAL.OFTInteger)
            #     # Pos attributes should not be required, the graph is noded (split at each intersection)
            # end

            objid = 0
            for restric in turn_restrictions
                # skip = false
                # for (fr, to) in zip(restric.geom[1:end - 1], restric.geom[2:end])
                #     if euclidean_distance(fr, to) < 1
                #         @error "Geometry for restriction $(restric.osm_id) has segments that are too short"
                #         skip = true
                #         break
                #     end
                # end
                # if skip  # and they say loop labels are an antipattern
                #     continue
                # end

                if length(restric.segments) ≤ 1
                    @error "Insufficient edges in $(restric.osm_id)"
                    continue
                end

                ArchGDAL.createfeature(layer) do f
                    lons = getproperty.(restric.geom, [:lon])
                    lats = getproperty.(restric.geom, [:lat])

                    ArchGDAL.setgeom!(f, ArchGDAL.createlinestring(lons, lats))
                    ArchGDAL.setfield!(f, 0, (objid += 1))  # ObjectID
                    ArchGDAL.setfield!(f, 1, restric.turn_angle) # Bearing
                    ArchGDAL.setfield!(f, 2, 1) # restricted
                    ArchGDAL.setfield!(f, 3, restric.osm_id) # OSM ID
                    # ArchGDAL.setfield!(f, 4, restric.edge_1_end ? 1 : 0)
                    # for (i, edge) in enumerate(restric.edges)
                    #     next_field = (i - 1) * 2 + 5
                    #     ArchGDAL.setfield!(f, next_field, edge)
                    #     ArchGDAL.setfield!(f, next_field + 1, 4)  # placeholder
                    # end
                end
            end

            ArchGDAL.copy(layer, dataset=ds)
        end
    end
end

"Find the way segments that precede and follow this node"
function find_way_segment_by_node(candidates::Vector{EdgeRef}, node_id::Int64)
    before_ref = nothing
    after_ref = nothing

    for candidate in candidates
        if candidate.nodes[end] == node_id
            isnothing(before_ref) || error("node $node_id appears more than once")
            before_ref = candidate
        end

        if candidate.nodes[1] == node_id
            isnothing(after_ref) || error("node $node_id appears more than once")
            after_ref = candidate
        end
    end

    # okay for one of them to be nothing
    isnothing(before_ref) && isnothing(after_ref) && error("node $node_id not found in way")

    return before_ref, after_ref
end

"Find a vertex in the center of a line segment"
function find_center(geom::AbstractVector{LatLon{Float64}}; frac=0.5)
    distances = Vector{Float64}()
    sizehint!(distances, length(geom))
    push!(distances, 0)
    for (fr, to) in zip(geom[1:end - 1], geom[2:end])
        seg_len = euclidean_distance(fr, to)
        push!(distances, distances[end] + seg_len)
    end

    # find the first vertex after the middle
    target_dist = distances[end] * frac
    first_after = findfirst(distances .≥ target_dist)
    last_before = first_after - 1

    # find the fraction along this line segment
    seg_frac = (target_dist - distances[last_before]) / (distances[first_after] - distances[last_before])
    g1 = geom[last_before]
    g2 = geom[first_after]
    return LatLon(g1.lat + (g2.lat - g1.lat) * seg_frac, g1.lon + (g2.lon - g1.lon) * seg_frac)
end


"""
Create an ArcGIS-format geometry for a turn restriction.

From https://desktop.arcgis.com/en/arcmap/latest/extensions/network-analyst/creating-a-turn-feature.htm

"Sequentially click each line feature that makes up the turn. Multiple vertices can be placed along a single edge element,
but at least one vertex must be placed on each edge element of the turn."

For U-turns on a single edge, things are a bit more complicated. From
https://desktop.arcgis.com/en/arcmap/latest/extensions/network-analyst/creating-a-turn-feature-representing-u-turns.htm

"Click somewhere along the length of the line feature.
Click at the end of the line feature where the U-turn occurs.
Double-click along the length of the same line feature to finish the sketch."

So we need three points with one at the vertex, contrary to usual instructions where you can't
have turn features at intersections (not documented but Arc will give an error when importing.)
"""
function create_turn_geometry(segments::AbstractVector{EdgeRef}, back::AbstractVector{Bool}, node_locations::Dict{Int64, LatLon})::AbstractVector{LatLon}
    if length(segments) != 2 || segments[1].way != segments[2].way
        # place two vertices on each edge, so that turn direction is clear in a feature
        # such as https://www.openstreetmap.org/relation/3644873
        geom = map(zip(segments, back)) do (seg, segback)
            geom = get.([node_locations], seg.nodes, [nothing])
            return [
                find_center(geom; frac=segback ? 0.6 : 0.4),
                find_center(geom; frac=segback ? 0.4 : 0.6)
            ]
        end |>
        Iterators.flatten |>
        collect
    else
        # special case for U-turn
        geom = [
            find_center(get.([node_locations], segments[1].nodes, [nothing])),
            node_locations[segments[1].nodes[back[1] ? 1 : end]],
            find_center(get.([node_locations], segments[2].nodes, [nothing]))
        ]
    end

    return geom
end

function process_simple_restriction(restric, from, to, via, way_segment_idx, node_locations)
    # all the segments each way got split into
    if !haskey(way_segment_idx, from.id)
        @warn "Way $(from.id) referenced in restriction $(restric.id), but not in graph"
        return
    end

    if !haskey(way_segment_idx, to.id)
        @warn "Way $(to.id) referenced in restriction $(restric.id), but not in graph"
        return
    end

    origin_candidates = way_segment_idx[from.id]
    destination_candidates = way_segment_idx[to.id]
    
    # locate the correct way segment
    local from_before, from_after, to_before, to_after
    try
        from_before, from_after = find_way_segment_by_node(origin_candidates, via.id)
    catch e
        @warn "Processing from way $(from.id) for restriction $(restric.id), $(e), skipping this restriction. Stacktrace:\n" *
            join(string.(stacktrace(catch_backtrace())), "\n")
        return nothing
    end
    
    try
        to_before, to_after = find_way_segment_by_node(destination_candidates, via.id)
    catch e
        @warn "Processing from way $(from.id) for restriction $(restric.id), $(e), skipping this restriction. Stacktrace:\n" *
            join(string.(stacktrace(catch_backtrace())), "\n")
        return nothing
    end
    
    #= check all possible combinations
    this is complicated by the fact that OSM ways are undirected. Consider a situation like this:
    
    ---------Main St--------+-------------------------
                            |
                            E                ^
                            l                N
                            m

                            S
                            t
                            |
                            |

    Suppose there is a no-left restriction from Main St to Elm St. We have to actually look at
    turn angles to figure out that this means WB Main -> SB Elm is prohibited, and EB Main -> SB Elm
    is allowed.

    It might be useful to have additional tags in restrictions for from_direction = forward/backward and to_direction
    # forward on both edges
    
    =#

    rtype = get_rtype(restric)

    # if both ways begin/end at the via node, then it is clear what was (probably) intended
    # ⊻ is exclusive or
    if (isnothing(from_before) ⊻ isnothing(from_after)) && (isnothing(to_before) ⊻ isnothing(to_after))
        from = isnothing(from_before) ? from_after : from_before
        to = isnothing(to_before) ? to_after : to_before

        # we are approaching the turn from the part of the way after the node, so traversing backwards
        from_back = isnothing(from_before)

        # we are leaving the turn from the part of the way before the node, do traversing backwards
        to_back = isnothing(to_after)

        bearing = get_turn_angle(from, from_back, to, to_back, node_locations)

        try
            if !is_turn_type(from, from_back, to, to_back, rtype, node_locations)
                @warn "Restriction $(restric.id) is nonambiguous, but indicates it should be of type $(rtype), but has bearing $(bearing)°"
            end
        catch e
            @warn "Could not process $(restric.id) turn direction, but restriction is nonambiguous: " *
                string(e) *
                join(stacktrace(catch_backtrace()), "\n")
        end

        return TurnRestriction(
            [from, to],
            [from_back, to_back],
            create_turn_geometry([from, to], [from_back, to_back], node_locations),
            bearing,
            restric.id
        )
    else
        candidates = [
            (from_before, false, to_after, false),  # traverse both forwards
            (from_before, false, to_before, true),
            (from_after, true, to_after, false),
            (from_after, true, to_before, true)
        ]

        remaining_candidates = filter(candidates) do c
            !isnothing(c[1]) && !isnothing(c[3]) && is_turn_type(c..., rtype, node_locations)
        end

        if length(remaining_candidates) == 0
            @warn "Processing restriction $(restric.id), no possible turns correspond to a $rtype direction, dropping this restriction. Possible turns:\n" *
            join(map(candidates) do c
                if !isnothing(c[1]) && !isnothing(c[2])
                    bearing = get_turn_angle(c..., node_locations)
                    return " - way $(from.id) $(c[2] ? backwards : forwards) => way $(to.id) $(c[4] ? backwards : forwards) " +
                    "($(bearing)°)"
                else
                    return nothing
                end
            end, "\n")
            return
        elseif length(remaining_candidates) > 1
            backwards = "backwards"
            forwards = "forwards"
            @warn "Restriction $(restric.id) is ambiguous, skipping. Possible candidates:\n" *
                join(map(remaining_candidates) do c
                    bearing = get_turn_angle(c...)
                    return " - way $(from.id) $(c[2] ? backwards : forwards) => way $(to.id) $(c[4] ? backwards : forwards) " +
                    "($(bearing)°)"
                end, "\n")
            return nothing
        else
            @info "generating restriction for $(restric.id)"
            from, from_back, to, to_back = remaining_candidates[1]
            bearing = get_turn_angle(from, from_back, to, to_back, node_locations)
    
            return TurnRestriction(
                [from, to],
                [from_back, to_back],
                create_turn_geometry([from, to], [from_back, to_back], node_locations),
                bearing,
                restric.id
            )
        end
    end
end

"""
Process a complex turn restriction with via way(s)
"""
function process_complex_restriction(restric, from, to, way_segment_idx, node_locations)
    # I don't think that the members of a relation are ordered... so the order the
    # via ways occur may or may not match the order they are traversed 🤦. Do a little depth-first search
    # to find all the ways you could hook them up.

    # because the way segments rather than the nodes are what we care about, we use a "dual graph"
    # with nodes representing way segments and way segments representing nodes. See
    # Winter, Stephan. 2002. “Modeling Costs of Turns in Route Planning.” GeoInformatica 6
    # (4): 345–361. doi:10.1023/A:1020853410145.

    edges_for_node = Dict{Int64, Vector{EdgeRef}}()
    wayid_for_edgeid = Dict{EdgeRef, Int64}()
    edges = Vector{EdgeRef}()
    ways = Vector{Int64}()

    for member in restric.members
        if !(member.role in Set(["from", "to", "via"]))
            continue  # warning printed above
        end

        if !haskey(way_segment_idx, member.id)
            @warn "processing restriction $(restric.id), way $(member.id) not found in graph"
            return nothing
        end

        segs = way_segment_idx[member.id]
        for seg in segs
            # graph is fully noded, so we only need to index first and last node
            for node in [seg.nodes[1], seg.nodes[end]]
                if !haskey(edges_for_node, node)
                    edges_for_node[node] = Vector{Int32}()
                end
                push!(edges_for_node[node], seg)
            end

            wayid_for_edgeid[seg] = member.id
            push!(edges, seg)
        end
        push!(ways, member.id)
    end

    g = Graph(length(edges))
    vidx_for_edge = Dict{EdgeRef, Int64}()
    for (i, edge) in enumerate(edges)
        # because it is a dual graph edge segment properties get set on vertices
        vidx_for_edge[edge] = i
    end

    for (node, edges) in edges_for_node
        # add edges to the dual graph that represent the connections made by this node
        for (efrom, eto) in Iterators.product(edges, edges)
            if efrom != eto
                add_edge!(g, vidx_for_edge[efrom], vidx_for_edge[eto])
            end
        end
    end

    # perform the dfs
    # figure out the origins and destinations
    origins = Set(map(s -> vidx_for_edge[s], way_segment_idx[from.id]))
    dests = Set(map(s -> vidx_for_edge[s], way_segment_idx[to.id]))

    edge_for_vidx = Dict(map(p -> p.second=>p.first, collect(vidx_for_edge))...)

    all_paths::Vector{Vector{EdgeRef}} = map(find_paths(g, origins, dests)) do path
        return map(path) do vx
            edge_for_vidx[vx]
        end
    end
    
    # paths that pass through all mentioned ways
    candidate_paths = filter(all_paths) do path
        # figure out if it traverses all the member ways
        traverses = Set{Int64}()
        for vx in path
            push!(traverses, wayid_for_edgeid[vx])
        end

        # make sure it traversed all ways
        for way in ways
            if !(way in traverses)
                return false
            end
        end
        return true
    end

    rtype = get_rtype(restric)

    # if there's only one candidate path, then it's easy
    if length(candidate_paths) == 1
        parsed = construct_restriction_from_path(candidate_paths[1], restric.id, node_locations)
        if !is_turn_type(parsed.turn_angle, rtype)
            @warn "Restric $(restric.id) has type $rtype but has bearing $(parsed.turn_angle)"
        end
        return parsed
    else
        # multiple candidate paths; filter by turn type
        candidate_restrictions =
            filter(map(path -> construct_restriction_from_path(path, restric.id, node_locations), candidate_paths)) do restriction
                is_turn_type(restriction.turn_angle, rtype)
            end
        
        if length(candidate_restrictions) == 1
            return candidate_restrictions[1]
        elseif length(candidate_restrictions) > 1
            @error "Restriction $(restric.id) is ambiguous, multiple possible restricted paths:\n" *
                join(map(candidate_restrictions) do r
                    return " - " * join(map(s -> "$(s.way) [$(s.nodes[1]) to $(s.nodes[end])]", r.segments), " → ")
                end, "\n")
            return nothing
        elseif length(candidate_restrictions) == 0
            @error "No matching restriction for $rtype, for restriction $(restric.id)"
            return nothing
        end
    end
end

function construct_restriction_from_path(path::Vector{EdgeRef}, restric_id, node_locations)
    back = falses(length(path))
    for (i, j) in zip(1:(length(path) - 1), 2:length(path))
        ei = path[i]
        ej = path[j]

        # should be exactly one duplicated node where they connect
        # if not, cannot determine forward/back
        nnodes = length(Set([ei.nodes[1], ei.nodes[end], ej.nodes[1], ej.nodes[end]]))
        if nnodes != 3
            error("path directions ambiguous, expected 3 nodes but got $nnodes")
        end

        found = false
        if ei.nodes[end] == ej.nodes[1]
            back[i] = back[j] = false
            found = true
        end

        if ei.nodes[1] == ej.nodes[1]
            found && error("multiple path directions found")
            back[i] = true
            back[j] = false
            found = true
        end
        
        if ei.nodes[end] == ej.nodes[end]
            found && error("multiple path directions found")
            back[i] = false
            back[j] = true
            found = true
        end

        if ei.nodes[1] == ej.nodes[end]
            found && error("multiple path directions found")
            back[i] = back[j] = true
            found = true
        end

        !found && error("unreachable state")
    end

    bearing = get_turn_angle(path[1], back[1], path[end], back[end], node_locations)

    return TurnRestriction(
        path,
        back,  # if edge 1 is back, then it does not pass through the end
        create_turn_geometry(path, back, node_locations),
        bearing,
        restric_id
    )
end

function get_turn_angle(first, first_back, second, second_back, node_locations)
    first_in, first_out = calculate_headings_for_geom(get.([node_locations], first.nodes, nothing))
    second_in, second_out = calculate_headings_for_geom(get.([node_locations], second.nodes, nothing))

    in_hdg = first_back ? (first_in + 180) % 360 : first_out
    out_hdg = second_back ? (second_out + 180) % 360 : second_in

    return bearing_between(in_hdg, out_hdg)
end

# Is the angle between edges (traveresed in the forward direction if _back is false, otherwise reverse)
# consistent with the restriction type listed in type?
function is_turn_type(first, first_back, second, second_back, type, node_locations)
    bearing = get_turn_angle(first, first_back, second, second_back, node_locations)

    return is_turn_type(bearing, type)
end

function is_turn_type(bearing, type)
    if endswith(type, "left_turn")
        target_ranges = LEFT_TURN_RANGES
    elseif endswith(type, "right_turn")
        target_ranges = RIGHT_TURN_RANGES
    elseif endswith(type, "u_turn")
        target_ranges = U_TURN_RANGES
    elseif endswith(type, "straight_on")
        target_ranges = STRAIGHT_RANGES
    else
        error("Unrecognized turn restriction $type")
    end

    for range in target_ranges
        # ranges will overlap at ends b/c using ≤ and ≥, but that's okay
        if range[1] ≤ bearing && range[2] ≥ bearing
            return true
        end
    end
    return false # not in any target range
end
