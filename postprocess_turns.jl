#=
postprocess the turn restrictions so they comply with Arc's arbitrary
restrictions on turn restrictions: https://pro.arcgis.com/en/pro-app/2.7/help/analysis/networks/network-dataset-build-errors.htm
specifically, clean turns that
- are duplicated
- One of the interior edges of the turn element is the same as the first or last edge.
- The edges of the turn element conflict with existing interior/exterior edges - interior
  elements cannot be exterior elements of other turns.
- do not have edges that are self-loops (TODO need to do this in network build step)

Fixing these has varying degrees of difficulty. Duplicated turns are easy - just remove them (warn, but probably a data error).

Interior edges is the same as the first or last edge - this *should* only happen with the last edge in restrictions
generated from an only_* restriction, I think.                 There are situations where you want this. For instance, consider the network in the docstring of only_turn.jl - the maneuver 
EB Cypress -> NB Elm -> SB Elm via u-turn @ Oak has the same (bidirectional) segment of Elm both in the middle and at the end.

A more realistic real-world example is here: https://www.openstreetmap.org/relation/7644917 . The generated turn restrictions include
one that has the from and via way and then a u-turn back onto the via way. ArcGIS can't handle this.

The workaround we use here is to create a restriction on the U-turn itself, without the restriction starting at the from way.
This is not strictly correct, as it will disallow for example NB Elm to SB Elm via U-turn at Oak in the example in the docstring,
but in most practical situations (for instance the one linked above) will suffice - and worst-case we ban some sketchy
U-turns that happen in the middle of turn restrictions. TODO never mind - will use turn restriction expansion defined below.

Edges of the turn element conflict with existing interior/exterior edge. First we need to index all interior edges. Then, for any turn
that starts or ends with one of these interior edges, we replace it with turns from all possible predecessors or to all possible descendants
of the problematic edge, so the problematic edge becomes an interior edge). Hopefully this doesn't get too out of hand - I think this generally
comes up when you have a no-U-turn restriction in the middle of a complex only_ restriction, that is fixed as above. (NB this should be used to
address the U-turn problem above, rather than breaking U-turns)

Self-loops are easy to fix by just eliminating self-loops in the graph build process. Make sure to add no-U-turn restrictions where
we break self-loops.

Note that one of these cleaning steps may cause more problems in the graph, so need to keep running repeatedly until there are no errors.
=#

"""
Expand a turn restriction, either at the start or at the end, by creating turn restrictions for all possible maneuvers
that could precede or follow said turn restriction
"""
function expand_turn_restriction(restric::TurnRestriction, start::Bool,
  node_locations::Dict{Int64, LatLon}, edges_for_node::Dict{Int64, Vector{EdgeRef}})

  # first node is end of first edge if back
  first_node = restric.segments[1].nodes[restric.back[1] ? end : 1]
  # last node is first node of last edge if back
  last_node = restric.segments[end].nodes[restric.back[end] ? 1 : end]

  node_to_expand = start ? first_node : last_node
  edge_to_expand = start ? restric.segments[1] : restric.segments[end]

  edges_to_expand = filter(e -> e !== edge_to_expand, edges_for_node[node_to_expand])
  # should remove exactly one edge
  @assert length(edges_to_expand) == length(edges_for_node[node_to_expand]) - 1

  return map(edges_to_expand) do edge
    if start
      edge_back = edge.nodes[1] == first_node
      edges = [edge; restric.segments]
      back = [edge_back; restric.back]
    else
      edge_back = edge.nodes[end] == last_node
      edges = [restric.segments; edge]
      back = [restric.back; edge_back]
    end

    geom = create_turn_geometry(edges, back, node_locations)
    bearing = get_turn_angle(edges[1], back[1], edges[end], back[end], node_locations)

    return TurnRestriction(
      edges,
      back,
      geom,
      bearing,
      restric.osm_id
    )
  end
end

function index_interior_edges(turns::Vector{TurnRestriction})
  interior_edges = Set()
  for turn in turns
    for edge in turn.segments[2:end - 1]
      push!(interior_edges, edge)
    end
  end
  return interior_edges
end

function postprocess_turns(turns::Vector{TurnRestriction}, node_locations::Dict{Int64, LatLon}, edges_for_node::Dict{Int64, Vector{EdgeRef}})
  iter = 0
  while true
    iter += 1
    @info "Postprocess turn restrictions: begin iteration $iter"
    
    # create a set of all interior edges
    interior_edges = index_interior_edges(turns)
    new_turns = Vector{TurnRestriction}()
    sizehint!(new_turns, length(turns))

    fixed_turns = 0
    removed_turns = 0

    unique_turns = Set()

    for turn in turns
      # when checking for uniqueness, only check edge sequence and forward/back,
      # don't check osm_id as there may be duplicate turns in OSM
      turn_summary = (turn.segments, turn.back)
      if turn_summary in unique_turns
        removed_turns += 1
        continue
      else
        push!(unique_turns, turn_summary)
      end

      # check for turns that start in the interior of a turn (themself or others)
      if turn.segments[1] in interior_edges
        # this does affect interior edges, but we loop through all turns, then re-index interior edges
        # and iterate
        append!(new_turns, expand_turn_restriction(turn, true, node_locations, edges_for_node))
        fixed_turns += 1
      elseif turn.segments[end] in interior_edges
        # elseif is intentional. it is very possible that a turn would start and end in the
        # interior of a turn, but we have already pushed the front-expanded turns above. if it's
        # both, we'll catch it on the next iteration
        append!(new_turns, expand_turn_restriction(turn, false, node_locations, edges_for_node))
        fixed_turns += 1
      else
        # this turn is okay
        push!(new_turns, turn)
      end
    end

    if fixed_turns == 0 && removed_turns == 0
      @info "All turn meet ArcGIS standards after $iter iteration(s)"
      return new_turns
    else
      n_new = length(new_turns) + removed_turns - length(turns)
      @info "In iteration $iter, expanded $fixed_turns turns into $n_new turns, and removed $removed_turns duplicate turns"
      turns = new_turns
    end
  end
end