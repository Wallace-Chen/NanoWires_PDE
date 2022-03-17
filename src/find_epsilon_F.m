% calculate the Fermi level
function epsilon_F = find_epsilon_F( V_input, E_schrodinger, psi_sqrt_m )

    global triangle_areas_poiss;
    global p_poiss; global t_poiss;
    
    % evaluate values at triangle midpoints, a row vector
    V_midpoint = pdeintrp( p_poiss, t_poiss, V_input );
    
    epsilon_F = fzero(@(E_F)charge_neutral( E_F, V_midpoint, triangle_areas_poiss, E_schrodinger, psi_sqrt_m), E_schrodinger(end) );

end