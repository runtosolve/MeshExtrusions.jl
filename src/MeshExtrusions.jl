module MeshExtrusions

using CrossSectionGeometry, LinesCurvesNodes 

function open_cross_section_with_solid_elements(centerline_XY_coordinates::Vector{Vector{Float64}}, outside_XY_coordinates::Vector{Vector{Float64}}, inside_XY_coordinates::Vector{Vector{Float64}}, Z::Vector{Float64})

    #Define number of solid element layers through the thickness:  (only 2 available right now)
    num_solid_element_layers=2

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

    XY = [XY; stack(centerline_XY_coordinates, dims=1)]
    XY = [XY; stack(inside_XY_coordinates, dims=1)]

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
        nodes = vcat(nodes, XYZ)

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

    top_surface_elements = Int[]

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

                if j == 1
                    push!(top_surface_elements, element_number)
                end

                solid_element_cell = [element_number n_1 n_2 n_3 n_4 n_5 n_6 n_7 n_8]

                solid_elements[element_number, :] = solid_element_cell

            end

        end

    end

    return nodes, solid_elements, top_surface_elements

end






function open_cross_section_with_shell_elements(X, Y, Z)


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








function closed_cross_section_with_shell_elements(X, Y, Z)


    num_cross_sections = length(Z)
    num_nodes_cross_section  = length(X)
    cross_section_nodes = [X Y]
     
    nodes = Array{Float64}(undef, 0, 3)
    for i = 1:num_cross_sections

        cross_section_nodes_xyz = [cross_section_nodes ones(Float64, num_nodes_cross_section) * Z[i]]
        nodes = vcat(nodes, cross_section_nodes_xyz)

    end

    #Calculate the number of elements in the cross-section.
    num_elements_cross_section = num_nodes_cross_section

    #Define the number of model segments along the length.
    num_model_segments = num_cross_sections - 1

    #Calculate total number of elements in the model.
    num_elements = num_elements_cross_section * num_model_segments

    #Initialize the solid element definition array.

    shell_elements = zeros(Int64, (num_elements, 9))

    element_number = 0

    for j = 1:num_model_segments
        
        for i = 1:num_elements_cross_section

            if i != num_elements_cross_section

                n_1 = i + (j - 1) * num_nodes_cross_section
                n_2 = i + j * num_nodes_cross_section
                n_3 = n_2 + 1
                n_4 = n_1 + 1
                n_5 = 0
                n_6 = 0
                n_7 = 0
                n_8 = 0

            else 

                n_1 = i + (j - 1) * num_nodes_cross_section 
                n_2 = i + j * num_nodes_cross_section
                n_3 = 1 + j * num_nodes_cross_section
                n_4 = 1 + (j - 1) * num_nodes_cross_section
                n_5 = 0
                n_6 = 0
                n_7 = 0
                n_8 = 0

            end

            element_number = element_number + 1

            shell_element_cell = [element_number n_1 n_2 n_3 n_4 n_5 n_6 n_7 n_8]

            shell_elements[element_number, :] = shell_element_cell

        end

    end

    nodes = Union{Float64, Int64}[1:size(nodes)[1] nodes]

    return nodes, shell_elements

end




function closed_cross_section_with_solid_elements(centerline_XY_coordinates::Vector{Vector{Float64}}, outside_XY_coordinates::Vector{Vector{Float64}}, inside_XY_coordinates::Vector{Vector{Float64}}, Z::Vector{Float64})

    # num_solid_element_layers = 2
    # num_cross_sections = length(z)
    # num_cross_section_nodes_layer = size(cross_section, 1)
    # num_nodes_cross_section = (num_solid_element_layers + 1) * num_cross_section_nodes_layer
    # num_cross_section_nodes = num_nodes_cross_section

    num_solid_element_layers = 2

    #Define number of cross-sections:
    num_cross_sections = length(Z)

    #Define number of nodes in one node layer of cross-section:
    num_cross_section_nodes_layer = size(centerline_XY_coordinates, 1)

    #Define number of nodes in cross-section:
    num_nodes_cross_section = (num_solid_element_layers + 1) * num_cross_section_nodes_layer


    num_cross_section_nodes = num_nodes_cross_section


    # unit_node_normals = CrossSection.Geometry.calculate_cross_section_unit_node_normals(cross_section)

    # Δ = t/2
    # top_coords = CrossSection.Geometry.get_coords_along_node_normals(cross_section, unit_node_normals, Δ)

    # Δ = -t/2
    # bottom_coords = CrossSection.Geometry.get_coords_along_node_normals(cross_section, unit_node_normals, Δ)

    # X = [top_coords[i][1] for i in eachindex(cross_section)]
    # Y = [top_coords[i][2] for i in eachindex(cross_section)]

    # cross_section_nodes = [X Y]

    # X = [cross_section[i][1] for i in eachindex(cross_section)]
    # Y = [cross_section[i][2] for i in eachindex(cross_section)]

    # cross_section_nodes = [cross_section_nodes; [X Y] ]

    # X = [bottom_coords[i][1] for i in eachindex(cross_section)]
    # Y = [bottom_coords[i][2] for i in eachindex(cross_section)]
    # cross_section_nodes = [cross_section_nodes; [X Y] ]
     
    # #Generate the nodal coordinate array for the whole model.

    # nodes = Array{Float64}(undef, 0, 3)
    # for i = 1:num_cross_sections

    #     cross_section_nodes_xyz = [cross_section_nodes ones(Float64, num_nodes_cross_section) * z[i]]
    #     nodes = vcat(nodes, cross_section_nodes_xyz)

    # end

    #Define number of solid element layers through the thickness:  (only 2 available right now)
    # num_solid_element_layers=2

    # #Define number of cross-sections:
    # num_cross_sections = length(Z)

    # #Define number of nodes in one node layer of cross-section:
    # num_cross_section_nodes_layer = size(centerline_XY_coordinates, 1)

    # #Define number of nodes in cross-section:
    # num_nodes_cross_section = (num_solid_element_layers + 1) * num_cross_section_nodes_layer

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

    XY = [XY; stack(centerline_XY_coordinates, dims=1)]
    XY = [XY; stack(inside_XY_coordinates, dims=1)]

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
        nodes = vcat(nodes, XYZ)

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

                if i != num_solid_elements_cross_section_layer

                    n_1 = i + (j - 1) * num_cross_section_nodes_layer + (k - 1) * num_cross_section_nodes
                    n_2 = i + j * num_cross_section_nodes_layer + (k - 1) * num_cross_section_nodes
                    n_3 = n_2 + 1
                    n_4 = n_1 + 1
                    n_5 = n_1 + num_cross_section_nodes
                    n_6 = n_2 + num_cross_section_nodes
                    n_7 = n_3 + num_cross_section_nodes
                    n_8 = n_4 + num_cross_section_nodes

                elseif i == num_solid_elements_cross_section_layer #for element that closes tube

                    n_1 = i + (j - 1) * num_cross_section_nodes_layer + (k - 1) * num_cross_section_nodes
                    n_2 = i + j * num_cross_section_nodes_layer + (k - 1) * num_cross_section_nodes
                    n_3 = 1 + j * num_cross_section_nodes_layer + (k - 1) * num_cross_section_nodes
                    n_4 = 1 + (j - 1) * num_cross_section_nodes_layer + (k - 1) * num_cross_section_nodes
                    n_5 = n_1 + num_cross_section_nodes
                    n_6 = n_2 + num_cross_section_nodes
                    n_7 = n_3 + num_cross_section_nodes
                    n_8 = n_4 + num_cross_section_nodes

                end

                element_number = element_number + 1

                solid_element_cell = [element_number n_1 n_2 n_3 n_4 n_5 n_6 n_7 n_8]

                solid_elements[element_number, :] = solid_element_cell

            end

        end

    end

    return nodes, solid_elements

end



function mesh_w_section(bf, d, t, member_length)

    ### top flange 
    L = [bf]
    θ = [0.0]
    n = [2]

    cross_section = CrossSectionGeometry.generate_thin_walled(L, θ, n)

    X = [cross_section[i][1] for i in eachindex(cross_section)]
    Y = [cross_section[i][2] for i in eachindex(cross_section)] .+ (d - t)

    Z_n = Int(member_length)
    Z = collect(range(0.0, member_length, Z_n + 1))

    nodes_tf, elements_tf = open_cross_section_with_shell_elements(X, Y, Z)


    #bottom flange 
    L = [bf]
    θ = [0.0]
    n = [2]

    cross_section = CrossSectionGeometry.generate_thin_walled(L, θ, n)

    X = [cross_section[i][1] for i in eachindex(cross_section)]
    Y = [cross_section[i][2] for i in eachindex(cross_section)]

    Z_n = Int(member_length)
    Z = collect(range(0.0, member_length, Z_n + 1))

    nodes_bf, elements_bf = open_cross_section_with_shell_elements(X, Y, Z)

    nodes_bf[:, 1] .+= 10000
    elements_bf[:, 1:5] .+= 10000



    #web
    L = [d - t]
    θ = [π/2]
    n = [1]

    cross_section = CrossSectionGeometry.generate_thin_walled(L, θ, n)

    X = [cross_section[i][1] for i in eachindex(cross_section)] .+ bf / 2 
    Y = [cross_section[i][2] for i in eachindex(cross_section)] 

    Z_n = Int(member_length)
    Z = collect(range(0.0, member_length, Z_n + 1))

    nodes_web, elements_web = open_cross_section_with_shell_elements(X, Y, Z)

    nodes_web[:, 1] .+= 20000
    elements_web[:, 1:5] .+= 20000


    nodes_all = [nodes_tf;  nodes_bf; nodes_web]
    elements_all = [elements_tf;  elements_bf; elements_web]

    # using CairoMakie 
    # scatter(nodes_frame_long[:, 2:4], markersize = 1.0)


    # index = findall(z_node -> z_node == 0.0, nodes_frame_long[:, 4])
    # scatter(nodes_frame_long[index, 2:3])



    #find nodes at centerline of top flange 

    nodes = nodes_tf

    xloc = bf / 2
    yloc = d - t 
    zloc = 0.0
    atol_x = 0.01
    atol_y = 0.01
    atol_z = 999999999.0 
    index = LinesCurvesNodes.find_nodes(nodes[:, 2:4], xloc, yloc, zloc, atol_x, atol_y, atol_z)
    node_set_top_flange = nodes[index, 1]  


    nodes = nodes_bf

    xloc = bf / 2
    yloc = 0.0
    zloc = 0.0
    atol_x = 0.01
    atol_y = 0.01
    atol_z = 999999999.0 
    index = LinesCurvesNodes.find_nodes(nodes[:, 2:4], xloc, yloc, zloc, atol_x, atol_y, atol_z)
    node_set_bottom_flange = nodes[index, 1]  


    nodes = nodes_web

    xloc = bf / 2
    yloc = 0.0
    zloc = 0.0
    atol_x = 0.01
    atol_y = 0.01
    atol_z = 999999999.0 
    index = LinesCurvesNodes.find_nodes(nodes[:, 2:4], xloc, yloc, zloc, atol_x, atol_y, atol_z)
    node_set_web_bottom = nodes[index, 1]  


    xloc = bf / 2
    yloc = d - t
    zloc = 0.0
    atol_x = 0.01
    atol_y = 0.01
    atol_z = 999999999.0 
    index = LinesCurvesNodes.find_nodes(nodes[:, 2:4], xloc, yloc, zloc, atol_x, atol_y, atol_z)
    node_set_web_top = nodes[index, 1]  



    node_set_remove = [node_set_web_bottom; node_set_web_top]

    for i in eachindex(node_set_remove)

        index = findfirst(node_number -> node_number == node_set_remove[i], nodes_all[:, 1])

        nodes_all = nodes_all[1:end .!= index, :]

    end

    element_nodes_all = elements_all[:, 2:end]

    for i in eachindex(node_set_web_bottom)

        indices = findall(node_number -> node_number == node_set_web_bottom[i], element_nodes_all)

        element_nodes_all[indices] .= node_set_bottom_flange[i]

    end


    for i in eachindex(node_set_web_top)

        indices = findall(node_number -> node_number == node_set_web_top[i], element_nodes_all)

        element_nodes_all[indices] .= node_set_top_flange[i]

    end


    elements_all = [elements_all[:, 1] element_nodes_all]


    # frame_y_nodes = nodes_all 
    # frame_y_elements = elements_all 

    return nodes_all, elements_all

end



end # module MeshExtrusions
