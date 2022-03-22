% solve the scalar schrodinger equation to retrieve at least N solutions
% a helper function for the iterate_run_new
function [psi_schrodinger, E_schrodinger, r_schrod_upper] = solve_eig(N, V_poiss, r_schrod_lower)

    global p_schrod; global e_schrod; global t_schrod;
    
    %% setup coefficient matrix
    % boundary matrix
    bmatrix = generate_dirichlet_boundary( 1 );
    c_schrod = 'c_coefficient_schrod(sd)';
    a_schrod = sprintf('V_total_without_xc(x, y, %s)', mat2str(V_poiss));
    d_schrod = 1;
    r_schrod = [r_schrod_lower, r_schrod_lower+1];
    
    psi_schrodinger = [];
    E_schrodinger = [];
    
    while length(E_schrodinger) < N
        [~, psi_schrodinger_add, E_schrodinger_add] = ...
            evalc('pdeeig(bmatrix, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod)');
        if ~isempty(E_schrodinger_add)
           psi_schrodinger = [psi_schrodinger, psi_schrodinger_add];
           E_schrodinger = [E_schrodinger; E_schrodinger_add];
        end
        r_schrod(1) = r_schrod(2);
        r_schrod(2) = r_schrod(2) + 1;
    end
    r_schrod_upper = r_schrod(1);
end