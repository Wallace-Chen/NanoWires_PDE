% function to calculate the 
function n_D_val = n_D_func(x, y, heaviside_n_D)

    global debug;
    global p_poiss;
    global t_poiss;
    global l0;
    global e_cgs;
    global eV;
    
    % an array-smart multiplication of heaviside_n_D and the user-supplied
    % doping density function evaluated at the finite element grid points; this
    % array-smart multiplication effectively produces a doping density function
    % that only contributes in spatial regions of the potential that are above
    % the Fermi level
    n_D_function = heaviside_n_D .* set_doping_density_function( p_poiss(1,:),p_poiss(2,:) );
    
    factor = 4 * pi * l0^2 * 1e-3 * e_cgs^2 / eV;
    n_D_val = factor * griddata_NaN( p_poiss(1,:), p_poiss(2,:), n_D_function, x, y );
    % n_D_val = factor * interpolant_NaN( p_poiss, t_poiss, n_D_function', x, y )';

end