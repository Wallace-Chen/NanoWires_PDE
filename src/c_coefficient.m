% a function to return c coefficient matrix for the PDE equation, c is
% constant for all equations, thus a scalar
function c = c_coefficient()

    global debug;

    global vector_of_masses;
    global vector_of_side_lengths;
    %%
    number_of_interfaces = length(vector_of_masses);
    
    c = [ mat2str(1 / vector_of_masses(1)) ];
    
    for i = 2:number_of_interfaces
        c = [c, ' + ', mat2str(1/vector_of_masses(i) - 1/vector_of_masses(i-1)), sprintf(' * heaviside_core_schrod(%f, x, y)', vector_of_side_lengths(i))];
    end
    
    %% debug purpose, showing the string representing the c matrix
    if debug
       disp("Debug: c_coefficient, showing the c coefficient matrix...");
       c
    end
    
end