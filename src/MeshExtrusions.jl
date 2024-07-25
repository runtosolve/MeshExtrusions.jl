module MeshExtrusions

using CrossSectionGeometry

function extrude_open_cross_section_with_solid_elements(centerline_XY_coordinates::Vector{Vector{Float64}}, outside_XY_coordinates::Vector{Vector{Float64}}, inside_XY_coordinates::Vector{Vector{Float64}}, Z::Vector{Float64}, num_solid_element_layers=2)

    #Define number of cross-sections:
    num_cross_sections = length(Z)

    #Define number of nodes in one node layer of cross-section:
    num_cross_section_nodes_layer = size(centerline_XY_coordinates, 1)

    #Define number of nodes in cross-section:
    num_nodes_cross_section = (num_solid_element_layers + 1) * num_cross_section_nodes_layer

    # #Find 
    # unit_node_normals = CrossSection.Geometry.calculate_cross_section_unit_node_normals(centerline_XY_coordinates)

    # Δ = member_thickness/2
    # top_coords = CrossSection.Geometry.get_coords_along_node_normals(centerline_XY_coordinates, unit_node_normals, Δ)

    # Δ = -member_thickness/2
    # bottom_coords = CrossSection.Geometry.get_coords_along_node_normals(centerline_XY_coordinates, unit_node_normals, Δ)

    XY = stack(outside_XY_coordinates, dims=1)

    # X = [top_coords[i][1] for i in eachindex(centerline_XY_coordinates)]
    # Y = [top_coords[i][2] for i in eachindex(centerline_XY_coordinates)]

    # cross_section_nodes = [X Y]

    XY = [XY; centerline_XY_coordinates]
    XY = [XY; inside_XY_coordinates]

    # X = [centerline_XY_coordinates[i][1] for i in eachindex(centerline_XY_coordinates)]
    # Y = [centerline_XY_coordinates[i][2] for i in eachindex(centerline_XY_coordinates)]

    # cross_section_nodes = [cross_section_nodes; [X Y] ]

    # X = [bottom_coords[i][1] for i in eachindex(centerline_XY_coordinates)]
    # Y = [bottom_coords[i][2] for i in eachindex(centerline_XY_coordinates)]
    # cross_section_nodes = [cross_section_nodes; [X Y] ]
    
    #Generate the nodal coordinate array for the whole model.

    nodes = Array{Float64}(undef, 0, 3)
    for i = 1:num_cross_sections

        # cross_section_nodes_xyz = [cross_section_nodes ones(Float64, num_nodes_cross_section) * Z[i]]
        XYZ = [XY ones(Float64, num_nodes_cross_section) * Z[i]]
        nodes = vcat(nodes, cross_section_nodes_xyz)

    end

    nodes = Union{Float64, Int64}[1:size(nodes)[1] nodes]

    #Calculate the number of solid elements in one layer of the cross-section.
    num_solid_elements_cross_section_layer = num_cross_section_nodes_layer - 1

    #Calculate the number of solid elements in the cross-section.
    num_solid_elements_cross_section = num_solid_elements_cross_section_layer * num_solid_element_layers

    #Define the number of model segments along the length.
    num_model_segments = num_cross_sections - 1

    #Calculate total number of solid elements in the model.
    num_solid_elements = num_solid_elements_cross_section * num_model_segments

    #Initialize the solid element definition array.

    solid_elements = zeros(Int64, (num_solid_elements, 9))

    element_number = 0

    for k = 1:num_model_segments

        for j = 1:num_solid_element_layers

            for i=1:num_solid_elements_cross_section_layer

                    n_1 = i + (j - 1) * num_cross_section_nodes_layer + (k - 1) * num_nodes_cross_section
                    n_2 = i + j * num_cross_section_nodes_layer + (k - 1) * num_nodes_cross_section
                    n_3 = n_2 + 1
                    n_4 = n_1 + 1
                    n_5 = n_1 + num_nodes_cross_section
                    n_6 = n_2 + num_nodes_cross_section
                    n_7 = n_3 + num_nodes_cross_section
                    n_8 = n_4 + num_nodes_cross_section

                element_number = element_number + 1

                solid_element_cell = [element_number n_1 n_2 n_3 n_4 n_5 n_6 n_7 n_8]

                solid_elements[element_number, :] = solid_element_cell

            end

        end

    end

end






function extrude_open_cross_section_with_shell_elements(X, Y, Z)


    num_cross_sections = length(Z)
    num_nodes_cross_section  = length(X)
    cross_section_nodes = [X Y]
     
    nodes = Array{Float64}(undef, 0, 3)
    for i = 1:num_cross_sections

        cross_section_nodes_xyz = [cross_section_nodes ones(Float64, num_nodes_cross_section) * Z[i]]
        nodes = vcat(nodes, cross_section_nodes_xyz)

    end

    #Calculate the number of elements in the cross-section.
    num_elements_cross_section = num_nodes_cross_section - 1

    #Define the number of model segments along the length.
    num_model_segments = num_cross_sections - 1

    #Calculate total number of elements in the model.
    num_elements = num_elements_cross_section * num_model_segments

    #Initialize the solid element definition array.

    shell_elements = zeros(Int64, (num_elements, 9))

    element_number = 0

    for j = 1:num_model_segments
        
        for i = 1:num_elements_cross_section

            n_1 = i + (j - 1) * num_nodes_cross_section
            n_2 = i + j * num_nodes_cross_section
            n_3 = n_2 + 1
            n_4 = n_1 + 1
            n_5 = 0
            n_6 = 0
            n_7 = 0
            n_8 = 0

            element_number = element_number + 1

            shell_element_cell = [element_number n_1 n_2 n_3 n_4 n_5 n_6 n_7 n_8]

            shell_elements[element_number, :] = shell_element_cell

        end

    end

    nodes = Union{Float64, Int64}[1:size(nodes)[1] nodes]

    return nodes, shell_elements

end


#for cross-sections with varying thickness, find the array range for each line segment that has the same thickness, used in CUFSM.Show
function find_linesegments(linewidths)

    index_start = 1

    linesegment_ranges = Vector{Vector{Int}}(undef, 0)
    linewidth_segments = Vector{Float64}(undef, 0)
    for i in eachindex(linewidths[1:end-1])

        if linewidths[i] != linewidths[i+1]

            index_end = i
            push!(linesegment_ranges, [index_start, index_end])
            push!(linewidth_segments, linewidths[i])
            index_start = i+1

        end

        if i == (length(linewidths) - 1)
            index_end = i + 1
            push!(linesegment_ranges, [index_start, index_end])
            push!(linewidth_segments, linewidths[i+1])
        end

    end

    return linesegment_ranges, linewidth_segments

end

#for cross-sections with varying thickness, used in CUFSM.Show
function combine_points_into_linesegments(linesegment_ranges, x, y)

    linesegments = Vector{Vector{Vector{Float64}}}(undef, size(linesegment_ranges)[1])
    # linesegments = Vector{Vector{Vector{Float64}}}(undef, size(linesegment_ranges)[1])

    for i in eachindex(linesegment_ranges)

        # linesegments[i] = cross_section_coords[linesegment_ranges[i][1]:linesegment_ranges[i][2]]
      
        nodes = Vector{Float64}[]
        for j = linesegment_ranges[i][1]:linesegment_ranges[i][2]

            push!(nodes, [x[j], y[j]])

            if (j+1) < length(x)  #stops 1 element short for a closed section, taken care of in Show.section to plot last closed section element
                push!(nodes, [x[j+1], y[j+1]])
            end

            # if (j+1) < length(x)  #stops 1 element short for a closed section, taken care of in Show.section to plot last closed section element
            # push!(nodes, [x[j+1], y[j+1]])
            # end

        end

   
        if (i == size(linesegment_ranges)[1]) & (length(x) > size(linesegment_ranges)[1])   #get last node in an open cross-section, for section plots
            push!(nodes, [x[end], y[end]])
        end
        

        nodes = unique(nodes)

        linesegments[i] = nodes


    end

    return linesegments

end

end # module MeshExtrusions
