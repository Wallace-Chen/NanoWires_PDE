% calculate the initial "seed" potential without any band bending.
function V_conduction_band_val = V_conduction_band(x, y)

    global debug;
    global vector_of_side_lengths;
    global vector_of_V0;
    global number_of_interfaces;
    
    global save_file;
    
    %%
    V_expression = [ mat2str(vector_of_V0(1)) ];
    for j = 2:number_of_interfaces
        V_expression = [V_expression, ' + ', mat2str( vector_of_V0(j)-vector_of_V0(j-1) ), sprintf(' * heaviside_core_schrod(%f, x, y)', vector_of_side_lengths(j))];
    end
   
    V_conduction_band_val = eval(V_expression);
    
    save(save_file, '-append');
    
    %% debug purpose, showing the string representing the v_cb
    if debug
        disp("Debug: V_conduction_band, the string representing the v_cb is:");
        V_expression
    end

end