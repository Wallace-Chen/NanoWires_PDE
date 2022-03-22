% solve electrostat and exchange potential from poisson equations

function [V_poiss, V_xc] = solve_potential(N, epsilon_F, E_schrodinger, psi_schrodinger, heaviside_n_D, factor)

    %% global variable declaration
    
    global debug;
    global p_poiss; global e_poiss; global t_poiss;
    global C;
    global psi_sqrt_eps_xc_m;

    a_poiss = 0;
    b_poiss = generate_dirichlet_boundary( 1 );
    c_poiss = 'c_coefficient_poiss(sd)';
    %% debug purpose
    %c_poiss = [ mat2str(-vector_of_eps(1)) ];
    % for j = 2 : number_of_interfaces
    %     c_poiss = [c_poiss, ' + ', mat2str(vector_of_eps(j-1) - vector_of_eps(j)), ...
    %         sprintf(' * heaviside_core_poiss(%s, x, y)', mat2str(vector_of_side_lengths(j)))];
    % end
    %% debug end
    
    % find normalization and |psi|^2*sqrt(m) term integrated over (x,y), normalized.
    [normalization, ~] = normalize_and_sqrt_m_triangular( psi_schrodinger );
    % calculate the psi_sqrt_eps term needed to calculate the total
    % electron density
    psi_sqrt_eps_sum = psi_sqrt_eps_summation( N, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    
    %% SCALAR PDE EQUATION
    % f coefficient for the scalar PDE equation, solve the electrostatic
    % potentials, due to electrons, dopants, and polarizations. V_cb do not
    % need to be solved (assumed to remain unchanged through the iterations).
    f_poiss_1 = sprintf('n_D_func(x, y, %s) + n_e_func(x, y, %s)', mat2str(heaviside_n_D), mat2str(psi_sqrt_eps_sum));
    if debug
       disp("debug: initiate_system, showing the coefficients for the scalar poisson equation...");
       c_poiss
       a_poiss
       f_poiss_1
    end
    V_poiss = assempde(b_poiss, p_poiss, e_poiss, t_poiss, c_poiss, a_poiss, f_poiss_1);
    V_poiss = factor * V_poiss / C;
    V_poiss = V_poiss - min(V_poiss, [], 'all');
    
    %% N*N PDE EQUATIONS, for exchange potentials
    psi_sqrt_eps_xc_m = psi_sqrt_eps_xc( N, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    % f_poiss_ex = f_ex_coefficient(N);
    % b_poiss_ex = generate_dirichlet_boundary( N*N );
    % if debug
    %    disp("debug: initiate_system, showing the coefficients for the N*N poisson equations for exchange potentials...");
    %   c_poiss
    %   a_poiss
    %   f_poiss_ex
    % end
    % compute the exchange potentials, this takes a lot of time
    % V_xc = assempde(b_poiss_ex, p_poiss, e_poiss, t_poiss, c_poiss, a_poiss, f_poiss_ex);
    % V_xc = V_xc / C;
    
    Np = length(V_poiss);
    V_xc = zeros( N*N*Np, 1 );
    for i = 1 : N*N
       V_xc( (i-1)*Np+1 : i*Np ) =  assempde(b_poiss, p_poiss, e_poiss, t_poiss, c_poiss, a_poiss, sprintf('n_e_xc_func(x, y, %d)', i) );
    end
    V_xc = factor * V_xc / C;
    V_xc = V_xc - min(V_xc, [], 'all');

end