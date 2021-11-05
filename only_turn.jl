"""
Using a restriction built for a no-turn feature, convert it to an only-turn feature,
by creating a series of restrictions that restrict everything other than the allowed turn.

This is straightforward for via nodes. For via ways it is more complicated. Consider an intersection
like this:

                      |
                      H        /
                      o       St
                      l      /
                      l     m 
                      y    l
                          E
                      St /
                      | /
                      |/
                      |
                      E
                      l
                      m

                      St
                      |
                      |----Oak St-------o (dead end)
---Cypress St---------|
                      |
                      |
                      |

Suppose that there is an only_left_turn restriction from EB Cypress St to NB Holly St via NB Elm St.
In this case, all of the following actions need to be restricted

EB Cypress -> SB Elm
EB Cypress -> NB Elm -> EB Oak
EB Cypress -> NB Elm -> NB Elm
EB Cypress -> NB Elm -> SB Elm via u-turn @ Oak  (tricky)
EB Cypress -> NB Elm -> SB Elm via u-turn @ Holly

So what needs to happen is for each of the edges in the restriction, we need to create a turn restriction starting with the
from way, passing through all the edges up to that point, and then ending with all edges _other_ than the one that allows you to continue.

It is critical that each restriction start at the very beginning of the only restriction and include all edges up to where things went off
the rails. Consider a counterexample. Suppose that we implemented the above restriction as:
No right from NB Cypress to SB Elm
No right from NB Elm to EB Oak
No straight from NB Elm to NB Elm at Holly
No U turn NB Elm at Oak

This would restrict the maneuvers above, but would also prevent these legal maneuvers:
NB Elm -> EB Oak
NB Elm -> NB Elm at Holly
qed

Only restrictions with via ways rather than nodes are uncommon, but do exist (e.g. https://www.openstreetmap.org/relation/7644917). This
algorithm is general enough to support arbitrary number of edges (though things could get out of hand if there are many).

It's technically undefined what an only restriction with a via way means. From the OSM wiki:  Going to other ways from the via point is forbidden with this relation.
and via is defined as a node. One reading would be to say that an only-right-turn restriction from ways A - B - C would mean that after
traveling A-B you can only go to C, but A-D would be allowed. But I think what this code actually does is more likely to be what is intended. From Way A at the node where
it intersects way B, you can only go B-C.
"""
function convert_restriction_to_only_turn(restric::TurnRestriction, node_locations::Dict{Int64, LatLon}, edges_for_node::Dict{Int64, Vector{EdgeRef}})
    restrictions = Vector{TurnRestriction}()

    for edgeidx in 1:(length(restric.segments) - 1)
        edge = restric.segments[edgeidx]
        next_edge = restric.segments[edgeidx + 1]
        # if it's a back-traversed edge, the end node is the start of the edge
        end_node = edge.nodes[restric.back[edgeidx] ? 1 : end]
        for eref in edges_for_node[end_node]
            if eref !== next_edge  # identity equality okay, should not be more than one of these per edge
                # create a restriction traversing all of the edges in the only restriction up to this point,
                # and then this disallowed maneuver. See above for why necessary.
                # if the next edge doesn't start with the end node it is a back edge
                eref.nodes[1] != end_node && eref.nodes[end] != end_node && error("neighbor list corrupted, contains edge not linked to node")
                eref_back = eref.nodes[1] != end_node
                segments = [restric.segments[1:edgeidx]; eref]
                back = [restric.back[1:edgeidx]; eref_back]

                geom = create_turn_geometry(segments, back, node_locations)
                bearing = get_turn_angle(segments[1], back[1], segments[end], back[end], node_locations)

                new_restric = TurnRestriction(
                    segments,
                    back,
                    geom,
                    bearing,
                    restric.osm_id
                )

                push!(restrictions, new_restric)
            end
        end
    end

    return restrictions
end