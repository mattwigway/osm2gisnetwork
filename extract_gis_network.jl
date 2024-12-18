# Extract ArcGIS network
# Extract a shapefile ready to be imported as an ArcGIS network dataset

using OpenStreetMapPBF
using ArchGDAL
using Logging
using ArgParse
using ProgressMeter
using Geodesy

EdgeRef = @NamedTuple{way::Int64, nodes::Vector{Int64}, oid::Int32, oneway::String}

include("compute_heading.jl")
include("graph_algos.jl")
include("process_turn_restrictions.jl")
include("only_turn.jl")
include("postprocess_turns.jl")

argtable = ArgParseSettings()
@add_arg_table! argtable begin
    "infile"
        help = "Input OSM PBF files"
    "outfile"
        help = "Output file"
    "--driver"
        help = "GDAL driver to use to write outfile, default GeoPackage (https://gdal.org/drivers/vector/index.html for full list, some may not be available on your system; must support multiple layers per file)"
        default = "GPKG"
end

const MILES_TO_KILOMETERS = 1.609344
const KNOTS_TO_KMH = 1.852

# car.lua says the traffic light penalty is 2 deciseconds or 0.2 seconds, but that seems unreasonably
# short. Based on reading the code it seems likely that these are intended to be _decaseconds_ or
# tens of seconds, and the docs are simply wrong
# https://github.com/Project-OSRM/osrm-backend/issues/5989
# const TRAFFIC_LIGHT_PENALTY_SECS = 20.0
# const BASE_TURN_PENALTY = 7.5
# # makes left turns more costly than right - for drive-on-left countries, set to
# # 1/1.075
# const TURN_BIAS = 1.075
# const TURN_PENALTY = 7.5

const DEFAULT_SPEED = 50

# copied from OSRM 5.24.0 car.lua
const freeflow_speeds_kmh = Dict(
    "motorway" => 90,
    "motorway_link" => 45,
    "trunk" => 85,
    "trunk_link" => 40,
    "primary" => 65,
    "primary_link" => 30,
    "secondary" => 55,
    "secondary_link" => 25,
    "tertiary" => 40,
    "tertiary_link" => 20,
    "unclassified" => 25,
    "residential" => 25,
    "living_street" => 10,
    "service" => 15,

    # additional link types that are driveable
    "road" => DEFAULT_SPEED
)

const hierarchies = Dict(
    "motorway" => 1,
    "motorway_link" => 1,
    "trunk" => 2,
    "trunk_link" => 2,
    "primary" => 2,
    "primary_link" => 2,
    "secondary" => 2,
    "secondary_link" => 2,
    "tertiary" => 2,
    "tertiary_link" => 2,
    "unclassified" => 3,
    "residential" => 3,
    "living_street" => 3,
    "service" => 3,

    # additional link types that are driveable
    "road" => 3
)

const DEFAULT_HIERARCHY = 3

function way_is_driveable(way::Way)::Bool
    if !haskey(way.tags, "highway") || !haskey(freeflow_speeds_kmh, way.tags["highway"])
        return false  # footway etc
    end

    if haskey(way.tags, "area") && way.tags["area"] == "yes"
        return false  # driveable area, often coincident with other roads
    end
    
    driveable = true
    if haskey(way.tags, "motor_vehicle")
        mv = way.tags["motor_vehicle"]
        if mv == "no" || mv == "agricultural" || mv == "forestry" ||
            mv == "agricultural;forestry" || mv == "delivery" || mv == "customers"
            # leaving private and permissive as some people live on private streets
            driveable = false
        end
    end

    if haskey(way.tags, "motorcar")
        mv = way.tags["motorcar"]
        if mv == "no" || mv == "agricultural" || mv == "forestry" ||
            mv == "agricultural;forestry" || mv == "delivery" || mv == "customers"
            # leaving private and permissive as some people live on private streets
            driveable = false
        end
    end

    return driveable
end




function parse_max_speed(speed_text)::Union{Float64, Missing}
    try
        return parse(Float64, speed_text)
    catch
        # not a raw km/hr number
        mtch = match(r"([0-9]+)(?: ?)([a-zA-Z/]+)", speed_text)
        if isnothing(mtch)
            @warn "unable to parse speed limit $speed_text"
            return missing
        else
            speed_scalar = parse(Float64, mtch.captures[1])
            units = lowercase(mtch.captures[2])

            if (units == "kph" || units == "km/h" || units == "kmph")
                return speed_scalar
            elseif units == "mph"
                return speed_scalar * MILES_TO_KILOMETERS
            elseif units == "knots"
                return speed_scalar * KNOTS_TO_KMH
            else
                @warn "unknown speed unit $units"
                return missing
            end
        end
    end
end

function speed_for_way(w::Way)::Float64
    # calculate travel time
    if haskey(w.tags, "maxspeed")
        maxspeed = parse_max_speed(w.tags["maxspeed"])
        if !ismissing(maxspeed) && isfinite(maxspeed)
            return maxspeed
        end
    end

    # fall through if no max speed tag or malformed

    highway = w.tags["highway"]
    if haskey(freeflow_speeds_kmh, highway)
        return freeflow_speeds_kmh[highway]
    else
        @warn "No speed type for $highway"
        return DEFAULT_SPEED
    end
end

function oneway(w::Way)::String
    # determine if a way is one way
    oneway = "No"
    if haskey(w.tags, "oneway")
        owt = w.tags["oneway"]
        if (owt == "yes" || owt == "1" || owt == "true")
            oneway = "FT"
        elseif (owt == "-1" || owt == "reverse")
            oneway = "TF"
        end
    end
    return oneway
end

function main()
    args = parse_args(argtable)
    infile = args["infile"]
    outfile = args["outfile"]
    driver = args["driver"]

    nodes_to_retain = Set{Int64}()
    intersection_nodes = Set{Int64}()

    total_ways = 0

    @info "Reading $infile and writing to $outfile"
    @info "Pass 1: Reading ways to identify nodes and intersections"
    waysp = ProgressUnknown()
    scan_pbf(infile, ways=way -> begin
        if way_is_driveable(way)
            for node in way.nodes
                if in(node, nodes_to_retain)
                    # if we've already seen this note in another way, it's an intersection node
                    push!(intersection_nodes, node)
                else
                    # otherwise just accumulate for geometry
                    push!(nodes_to_retain, node)
                end
            end
            ProgressMeter.next!(waysp)
            total_ways += 1
        end
    end)
    ProgressMeter.finish!(waysp)

    location_for_nodeid = Dict{Int64, LatLon}()
    edges_for_node = Dict{Int64, Vector{EdgeRef}}()

    expected_nodes = length(nodes_to_retain)
    @info "Pass 2: Reading $expected_nodes nodes"
    nodesp = Progress(expected_nodes)
    scan_pbf(infile, nodes=node -> begin
        if in(node.id, nodes_to_retain)
            # we don't need to save the whole node with tags, just the location
            location_for_nodeid[node.id] = LatLon(node.lat, node.lon)
            ProgressMeter.next!(nodesp)
        end
    end)
    ProgressMeter.finish!(nodesp)

    read_nodes = length(location_for_nodeid)
    
    expected_nodes == read_nodes || 
        error("Expected to find $(expected_nodes) nodes but found $(read_nodes)")

    @info "Read $read_nodes nodes"

    @info "Constructing GDAL geometries"
    way_segment_index = Dict{Int64, Vector{EdgeRef}}()

    # https://discourse.julialang.org/t/how-to-create-a-new-shapefile-containing-a-few-points/43454/3
    ArchGDAL.create(outfile, driver = ArchGDAL.getdriver(driver)) do ds
        # EPSG 4326 - WGS 84 - coordinate reference system used by OpenStreetMap
        ArchGDAL.createlayer(name="network", geom=ArchGDAL.wkbLineString, dataset=ds, spatialref=ArchGDAL.importEPSG(4326)) do layer
            ArchGDAL.addfielddefn!(layer, "OBJECTID", ArchGDAL.OFTInteger)
            ArchGDAL.addfielddefn!(layer, "highway", ArchGDAL.OFTString)
            ArchGDAL.addfielddefn!(layer, "name", ArchGDAL.OFTString)
            #ArchGDAL.addfielddefn!(layer, "ref", ArchGDAL.OFTString)
            ArchGDAL.addfielddefn!(layer, "way_id", ArchGDAL.OFTInteger64)
            ArchGDAL.addfielddefn!(layer, "fr_node_id", ArchGDAL.OFTInteger64)
            ArchGDAL.addfielddefn!(layer, "to_node_id", ArchGDAL.OFTInteger64)
            ArchGDAL.addfielddefn!(layer, "Minutes", ArchGDAL.OFTReal)
            ArchGDAL.addfielddefn!(layer, "Oneway", ArchGDAL.OFTString)
            ArchGDAL.addfielddefn!(layer, "hierarchy", ArchGDAL.OFTInteger)
            lats = Vector{Float64}()
            lons = Vector{Float64}()
            oid::Int32 = one(Int32)

            fprog = Progress(total_ways)
            scan_pbf(infile, ways = way -> begin
                if way_is_driveable(way)
                    ProgressMeter.next!(fprog)
                    if length(way.nodes) < 2
                        return
                    end
                    empty!(lats)
                    empty!(lons)
                    edges_for_way = Vector{EdgeRef}()
                    nodes = Vector{Int64}()
                    push!(nodes, way.nodes[1])

                    # prepopulate with first node
                    fr_node_id = way.nodes[1]
                    node = location_for_nodeid[fr_node_id]
                    push!(lats, node.lat)
                    push!(lons, node.lon)
                    prev_loc = node
                    length_meters = 0

                    for (i, nodeid) in enumerate(way.nodes[2:end])
                        if nodeid == nodes[end]
                            # skip duplicated nodes
                            continue
                        end

                        push!(nodes, nodeid)
                        node = location_for_nodeid[nodeid]
                        push!(lats, node.lat)
                        push!(lons, node.lon)

                        length_meters += euclidean_distance(prev_loc, node)
                        prev_loc = node

                        # +1, +2 b/c we skipped the first node
                        # break at intersection, end, or if next node would form a loop (ArcGIS can't handle turn restrictions
                        # that start or end on loop edges)
                        # Note that this does not handle loop edges with no intervening nodes, but that doesn't really make
                        # sense anyways
                        if in(nodeid, intersection_nodes) || (i + 1) == length(way.nodes) || way.nodes[i + 2] == fr_node_id
                            # break the way here, create a feature
                            ArchGDAL.addfeature(layer) do f
                                ArchGDAL.setgeom!(f, ArchGDAL.createlinestring(lons, lats))
                                ArchGDAL.setfield!(f, 0, oid)
                                ArchGDAL.setfield!(f, 1, way.tags["highway"])
                                if haskey(way.tags, "name")
                                    ArchGDAL.setfield!(f, 2, way.tags["name"])
                                else
                                    ArchGDAL.setfield!(f, 2, "")
                                end
            
                                ArchGDAL.setfield!(f, 3, way.id)
                                ArchGDAL.setfield!(f, 4, fr_node_id)
                                ArchGDAL.setfield!(f, 5, nodeid)

                                maxspeed = speed_for_way(way)
                                # note that this does not account for turn/intersection costs (accounted separately)
                                travel_time_minutes = (length_meters / 1000) / maxspeed * 60
                                ArchGDAL.setfield!(f, 6, travel_time_minutes)
                                
                                # figure out oneway
                                ArchGDAL.setfield!(f, 7, oneway(way))

                                ArchGDAL.setfield!(f, 8, get(hierarchies, way.tags["highway"], DEFAULT_HIERARCHY))

                                edge = (way=way.id, nodes=nodes, oid=oid, oneway=oneway(way))
                                push!(edges_for_way, edge)

                                if !haskey(edges_for_node, fr_node_id)
                                    edges_for_node[fr_node_id] = [edge]
                                else
                                    push!(edges_for_node[fr_node_id], edge)
                                end

                                if !haskey(edges_for_node, nodeid)
                                    edges_for_node[nodeid] = [edge]
                                else
                                    push!(edges_for_node[nodeid], edge)
                                end

                                # createfeature uses setfeature! instead of addfeature!, so fid needs to be defined
                                ArchGDAL.setfid!(f, oid)
                                oid += 1
                            end
                            # prepare for next way segment
                            empty!(lats)
                            empty!(lons)
                            length_meters = 0
                            push!(lats, node.lat)
                            push!(lons, node.lon)
                            nodes = Vector{Int64}()
                            fr_node_id = nodeid
                            push!(nodes, fr_node_id)
                        end
                    end

                    way_segment_index[way.id] = edges_for_way
                end
            end)
            ProgressMeter.finish!(fprog)
        end

        @info "Writing turn features"
        write_turn_features(infile, ds, way_segment_index, location_for_nodeid, edges_for_node)
    end
end

main()