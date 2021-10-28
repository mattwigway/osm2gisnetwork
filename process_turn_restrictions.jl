# create turn features
# schema:
# https://desktop.arcgis.com/en/arcmap/latest/extensions/network-analyst/turns-in-the-network-dataset.htm

using OSMPBF
using Logging
using ProgressMeter

# a turn is a left turn if between these headings
const LEFT_TURN_RANGES = [(-150, -30)]
const STRAIGHT_RANGES = [(-30, 30)]
const RIGHT_TURN_RANGES = [(30, 150)]
const U_TURN_RANGES = [(-Inf32, -150), (150, Inf32)]

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
    @info "Parsing turn restrictions"
    rprog = ProgressUnknown()
    scan_pbf(infile, relations=r -> begin
        # first, check if this is a turn restriction
        if haskey(r.tags, "type") && r.tags["type"] == "restriction" && any(haskey.([r.tags], ["restriction", "restriction:motorcar"]))
            ProgressMeter.next!(rprog)

            rtype = get_rtype(r)

            if startswith(rtype, "only_")
                @warn "Ignoring $rtype restriction, only_* restrictions not yet supported"
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
                    @info "parsing restriction with via node"
                    parsed = process_simple_restriction(r, from, to, via[1], way_segment_idx, node_locations)
                elseif length(via) ≥ 1 && all(map(v -> v.type == OSMPBF.way, via))
                    @info "parsing restriction with via way(s)"
                elseif length(via) == 0
                    @warn "restriction $(r.id) has no via members, skipping"
                else
                    @warn "via members of restriction $(r.id) are invalid (multiple nodes, mixed nodes/ways, relation members), skipping"
                end
            end
        end
    end)
end

"Find the way segments that precede and follow this node"
function find_way_segment_by_node(candidates::Vector{EdgeRef}, node_id::Int64)
    before_ref = nothing
    after_ref = nothing

    for candidate in candidates
        if candidate.nodes[1] == node_id
            isnothing(before_ref) || error("node $node_id appears more than once")
            before_ref = candidate
        end

        if candidate.nodes[end] == node_id
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
    
    found_combination = nothing

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

        try
            if !is_turn_type(from, from_back, to, to_back, rtype, node_locations)
                bearing = get_turn_angle(from, from_back, to, to_back, node_locations)
                @warn "Restriction $(restric.id) is nonambiguous, but indicates it should be of type $(rtype), but has bearing $(bearing)°"
            end
        catch e
            @warn "Could not process $(restric.id) turn direction, but restriction is nonambiguous: " *
                string(e) *
                join(stacktrace(catch_backtrace()), "\n")
        end

        @info "generating restriction for $(restric.id)"
        found_combination = from => to
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
        elseif length(remaining_candidates) > 1
            backwards = "backwards"
            forwards = "forwards"
            @warn "Restriction $(restric.id) is ambiguous, skipping. Possible candidates:\n" *
                join(map(remaining_candidates) do c
                    bearing = get_turn_angle(c...)
                    return " - way $(from.id) $(c[2] ? backwards : forwards) => way $(to.id) $(c[4] ? backwards : forwards) " +
                    "($(bearing)°)"
                end, "\n")
        else
            @info "generating restriction for $(restric.id)"
        end
    end
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

    bearing = get_turn_angle(first, first_back, second, second_back, node_locations)
    for range in target_ranges
        # ranges will overlap at ends b/c using ≤ and ≥, but that's okay
        if range[1] ≤ bearing && range[2] ≥ bearing
            return true
        end
    end
    return false # not in any target range
end