% calculate the source term of exchange potential as a function of inputted
% x and y values.
function n_e_xc = n_e_xc_func(x, y, idx)

    global debug;
    global p_schrod;
    global t_schrod;
    global l0;
    global e_cgs;
    global eV;
    global psi_sqrt_eps_xc_m;
    %%
       
    [x_sym, y_sym] = symmetrize_coordinates(x, y);

    factor = 4 * pi * e_cgs^2 / l0 / eV;
    n_e_xc = factor / 2 * abs( n_e_prefactor(x, y) .* ...
        griddata_NaN( p_schrod(1, :), p_schrod(2, :), psi_sqrt_eps_xc_m(:, idx), x_sym, y_sym ) );
        % interpolant_NaN( p_schrod, t_schrod, psi_sqrt_eps_xc_m(:, idx), x_sym, y_sym )' );
    
end