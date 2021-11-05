# Convert OpenStreetMap networks to ArcGIS networks

This tool converts OSM networks to ArcGIS-networks, in a two step process.

Steps:

1. Download OSM for your area of interest, in PBF format, e.g. from [GeoFabrik](https://download.geofabrik.de)
2. Run the tool: `julia --project extract_arcgis_network.jl osm_file.pbf output.shp` to create two network shapefiles, one with the network and one with turn restrictions. Write down the maximum number of edges printed at the end.
3. In ArcGIS, create a file geodatabase to hold the network (or use an existing file geodatabase)
4. Within the geodatabase, create a feature dataset to hold your network. Make sure the spatial reference is the same as the network shapefiles produced by step 2.
5. Use the "Feature Class to Feature Class" tool to copy the shapefiles into the feature dataset.
6. Enable Network Analyst in Customize -> Extensions
7. Right-click on your feature dataset and choose New -> Network dataset
8. Select that the main network feature class will participate in the network, but _not_ the turns
9. Accept defaults for the rest. Connectivity at the ends is correct, as the OSM ways have been broken down into separate line features at each intersection.
10. Create an "auto" mode (and other modes if desired)
11. After the network dataset is created, accept the prompt to build the network dataset
12. Use the "Create Turn Feature Class" tool to create a new turn feature class in the feature dataset. Set the number of edges to the number recorded in step 2, the input network dataset to the network dataset from step 7, and the template to the feature class in your dataset that contains turn features.
13. Add the original feature class containing the network, the original feature class with the turns, the network dataset, and the turn feature class to the map.
14. On the editing toolbar, start editing
15. Right-click on the original layer containing the turns, choose Selection -> Select all
16. Copy the selected features
17. Select the newly-created turn feature class
18. Paste the selected features into this feature class
19. Edit the network dataset, under properties, set up a new field "restriction" of type restriction with default prohibited.
20. Under "Evaluators" for the new restriction property, for the turns layer, set the evaluator to constant, use restriction.
21. In the travel modes tab, check "restriction" as a restriction to use.
22. If desired, add a global turn cost default evaluator to the minutes field
23. If desired, add a hierarchy to the attribute table - the shapefile output has a hierachy field with 1 for highways, 2 for secondary roads, and 3 for residential. Even if you do not use hierarchical routing, the hierarchy can be used to define turn costs
20. Rebuild the network dataset.