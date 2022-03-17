% evaluate the charge-neutrality constraint equation
function residual = charge_neutral( E_F, V_midpoint, triangle_areas, E_schrodinger, psi_sqrt_m )

    global p_poiss; global t_poiss;
    
    %%
    % the number of states occupied by electrons
    n_quanta = sum( E_F > E_schrodinger );
    
    % find the depletion region
    triangle_indices_above_E_F = V_midpoint > E_F;
    
    % column vector of n_D evaluated at mesh points
    n_D_val_evaluated_at_grid_points = set_doping_density_function( p_poiss(1, :), p_poiss(2, :) )';
    
    % evaluate average n_D for each mesh triangles, a column vector
    avg_n_D_val = 1/3 * ...
        (n_D_val_evaluated_at_grid_points( t_poiss(1, :) ) + ...
         n_D_val_evaluated_at_grid_points( t_poiss(2, :) ) + ...
         n_D_val_evaluated_at_grid_points( t_poiss(3, :) ));
    
    % n_D(x,y) integrated over depletion region
    residual = (triangle_indices_above_E_F .* triangle_areas) * avg_n_D_val;
    
    % take care of the contribution from electrons
    for i = 1 : n_quanta
        residual = residual - 1 / pi * sqrt( E_F-E_schrodinger(i) ) * psi_sqrt_m(i);
    end

end