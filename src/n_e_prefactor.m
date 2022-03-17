% multiplicative prefactor for calculating the electron density n_e:
% sqrt(m) / pi
function prefactor = n_e_prefactor(x, y)

    global vector_of_masses;
    global vector_of_side_lengths;
    global number_of_interfaces;

    prefactor_expression = [ mat2str(vector_of_masses(1)) ];
    for j = 2 : number_of_interfaces
        prefactor_expression = [prefactor_expression, ' + ', mat2str(vector_of_masses(j)-vector_of_masses(j-1)), ...
            sprintf(' * heaviside_core_schrod(%s, x, y)', mat2str(vector_of_side_lengths(j)))];
    end

    prefactor = 1 / pi * sqrt( eval(prefactor_expression) );
    
end