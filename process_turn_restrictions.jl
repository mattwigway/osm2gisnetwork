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
    edges::AbstractVector{Int32}
    segments::AbstractVector{EdgeRef}
    edge_1_end::Bool
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
        node_locations::Dict{Int64, LatLon})

    n_only = 0
    turn_restrictions = Vector{TurnRestriction}()
    @info "Parsing turn restrictions"
    rprog = ProgressUnknown()
    scan_pbf(infile, relations=r -> begin
        # first, check if this is a turn restriction
        if haskey(r.tags, "type") && r.tags["type"] == "restriction" && any(haskey.([r.tags], ["restriction", "restriction:motorcar"]))
            ProgressMeter.next!(rprog)

            rtype = get_rtype(r)

            if startswith(rtype, "only_") || rtype == "no_exit" || rtype == "no_entry"
                n_only += 1
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

                if length(via) == 1 && via[1].type == OSMPBF.node
                    parsed = process_simple_restriction(r, from, to, via[1], way_segment_idx, node_locations)
                    if !isnothing(parsed)
                        push!(turn_restrictions, parsed)
                    end
                elseif length(via) ≥ 1 && all(map(v -> v.type == OSMPBF.way, via))
                    parsed = process_complex_restriction(r, from, to, way_segment_idx, node_locations)
                    if !isnothing(parsed)
                        push!(turn_restrictions, parsed)
                    end
                elseif length(via) == 0
                    @warn "restriction $(r.id) has no via members, skipping"
                else
                    @warn "via members of restriction $(r.id) are invalid (multiple nodes, mixed nodes/ways, relation members), skipping"
                end
            end
        end
    end)

    @info "created $(length(turn_restrictions)) turn restrictions"

    if n_only > 0
        @warn "Ignored $n_only no_entry, no_exit, or only_* restrictions"
    end

    # find the maximum number of edges in turn restrictions
    max_edges = max(map(r -> length(r.edges), turn_restrictions)...)

    @info "all turn restrictions have $max_edges edges or less"

    # https://discourse.julialang.org/t/how-to-create-a-new-shapefile-containing-a-few-points/43454/3
    ArchGDAL.create(outfile, driver = ArchGDAL.getdriver("ESRI Shapefile")) do ds
        # EPSG 4326 - WGS 84 - coordinate reference system used by OpenStreetMap
        ArchGDAL.createlayer(geom=ArchGDAL.wkbLineString, spatialref=ArchGDAL.importEPSG(4326)) do layer
            ArchGDAL.addfielddefn!(layer, "ObjectID", ArchGDAL.OFTInteger)
            ArchGDAL.addfielddefn!(layer, "Bearing", ArchGDAL.OFTReal)
            ArchGDAL.addfielddefn!(layer, "OSMID", ArchGDAL.OFTInteger64)
            ArchGDAL.addfielddefn!(layer, "Edge1End", ArchGDAL.OFTInteger)
            for i in 1:max_edges
                ArchGDAL.addfielddefn!(layer, "Edge$(i)FID", ArchGDAL.OFTInteger)
                # Pos attributes should not be required, the graph is noded (split at each intersection)
            end

            objid = 0
            for restric in turn_restrictions
                ArchGDAL.createfeature(layer) do f
                    lons = getproperty.(restric.geom, [:lon])
                    lats = getproperty.(restric.geom, [:lat])

                    ArchGDAL.setgeom!(f, ArchGDAL.createlinestring(lons, lats))
                    ArchGDAL.setfield!(f, 0, (objid += 1))
                    ArchGDAL.setfield!(f, 1, restric.turn_angle)
                    ArchGDAL.setfield!(f, 2, restric.osm_id)
                    ArchGDAL.setfield!(f, 3, restric.edge_1_end ? 1 : 0)
                    for (i, edge) in enumerate(restric.edges)
                        ArchGDAL.setfield!(f, i + 3, edge)
                    end
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

"""
Process a simple turn restriction, with a via node
"""
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

        from_geom = get.([node_locations], from.nodes, nothing)
        if from_back
            reverse!(from_geom)
        end

        to_geom = get.([node_locations], to.nodes, nothing)
        if to_back
            reverse!(to_geom)
        end

        geom = [from_geom..., to_geom[2:end]...]

        return TurnRestriction(
            [from.fcid, to.fcid],
            [from, to],
            !from_back,  # if the first edge is traversed backwards, we do not go through the end of it
            geom,
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

            from_geom = get.([node_locations], from.nodes, nothing)
            if from_back
                reverse!(from_geom)
            end
    
            to_geom = get.([node_locations], to.nodes, nothing)
            if to_back
                reverse!(to_geom)
            end
    
            geom = [from_geom..., to_geom[2:end]...]
    
            return TurnRestriction(
                [from.fcid, to.fcid],
                [from, to],
                !from_back,  # if the first edge is traversed backwards, we do not go through the end of it
                geom,
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
    edges = map(e -> e.fcid, path)
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

    # make the geom
    geom = Vector{LatLon}()
    
    for i in 1:length(path)
        edge = path[i]
        seg_geom = get.([node_locations], edge.nodes, nothing)
        if back[i]
            reverse!(seg_geom)
        end

        # don't duplicate vertices at intersections
        append!(geom, i == 1 ? seg_geom : seg_geom[2:end])
    end

    bearing = get_turn_angle(path[1], back[1], path[end], back[end], node_locations)

    return TurnRestriction(
        edges,
        path,
        !back[1],  # if edge 1 is back, then it does not pass through the end
        geom,
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

struct PathState
    back::Union{PathState, Nothing}
    at_vertex::Int64
end

"return all acyclic paths from origin nodes to destination nodes"
function find_paths(g, origins, destinations)
    q = Queue{PathState}()
    states_at_dest = Vector{PathState}()

    # enqueue all origins
    for origin in origins
        enqueue!(q, PathState(nothing, origin))
    end

    # while there is anything on the queue, pop it off and explore from it
    # like dijkstra but without pqueue or vertex labels
    while length(q) > 0
        from_state = dequeue!(q)

        prev_vertices = Set{Int64}()
        back_state = from_state  # not from_state.back so we get current vertex in prev_vertices
        while !isnothing(back_state)
            push!(prev_vertices, back_state.at_vertex)
            back_state = back_state.back
        end

        for vertex in Graphs.all_neighbors(g, from_state.at_vertex)
              # don't loop, and only traverse a single origin edge (i.e. don't traverse all parts of an origin way) 
            if !(vertex in prev_vertices) && !(vertex in origins)
                next_state = PathState(from_state, vertex)
                if vertex in destinations
                    # we found a path
                    push!(states_at_dest, next_state)
                else
                    # not there yet
                    enqueue!(q, next_state)
                end
            end
        end
    end

    # convert states to something more useful by back-traversing
    paths = Vector{Vector{Int32}}()
    sizehint!(paths, length(states_at_dest))
    for state in states_at_dest
        path = Vector{Int32}()
        back_state = state
        while !isnothing(back_state)
            push!(path, back_state.at_vertex)
            back_state = back_state.back
        end
        reverse!(path)
        push!(paths, path)
    end

    return paths
end