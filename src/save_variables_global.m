
function save_variables_global(variables, v_names)

    %% the following variables will be saved
    global vector_of_V0;
    global vector_of_masses;
    global vector_of_eps;
    global number_of_triangles;
    global outer_shell_side_length;
    
    global avg_side_length_schrod;
    global avg_side_length_poiss;
    global p_poiss; global e_poiss; global t_poiss;
    global p_schrod; global e_schrod; global t_schrod;
    global number_mesh_points_schrod;
    global number_mesh_points_poiss;
    global triangle_indices_of_each_region_schrod;
    global triangle_indices_of_each_region_poiss;
    global triangle_areas_schrod;
    global triangle_areas_poiss;
    global C;

    %%
    global degree_of_polygon;
    global vector_of_side_lengths;
    global number_of_interfaces;
    
    name_of_save_file = sprintf('%d_charge_neutrality', degree_of_polygon);
    for i = 1 : number_of_interfaces
        name_of_save_file = [name_of_save_file, '_', num2str(vector_of_side_lengths(i))];
    end
    save([name_of_save_file, '.mat']);
    save([name_of_save_file, '.mat'], v_names{:}, '-append');

    fprintf("Variables saved to %s.\n", [name_of_save_file, '.mat']);
end