% function to solve poisson and schrodinger equations iteratively until the
% potential gets converged.
% @param V_elestat: initial guess, contains only eletrcostatic potential.
function [psi_schrodinger, E_schrodinger, epsilon_F, V_elestat, average_change] = ...
    iterate_run(psi_schrodinger, E_schrodinger, epsilon_F, V_elestat)

    global debug;
    global p_poiss;
    global e_poiss;
    global t_poiss;
    global p_schrod;
    global e_schrod;
    global t_schrod;
    global vector_of_eps;
    global number_of_interfaces;
    global vector_of_side_lengths;
    global C;
    
    % make it global to have neat code
    global psi_sqrt_eps_xc_m;
    %%
    % coefficient preparation for the poisson equations
    % c coefficient for the poisson equation, a scalar for whatever
    % dimension
%     c_poiss = [ mat2str(-vector_of_eps(1)) ];
%     for j = 2 : number_of_interfaces
%         c_poiss = [c_poiss, ' + ', mat2str(vector_of_eps(j-1) - vector_of_eps(j)), ...
%             sprintf(' * heaviside_core_poiss(%s, x, y)', mat2str(vector_of_side_lengths(j)))];
%     end
    c_poiss = 'c_coefficient_poiss(sd)';
    % a coefficient for the poisson equation, constant, equal to 0
    a_poiss = 0;
    % boundary dirichlet matrix for the scalar PDE equation
    b_poiss = generate_dirichlet_boundary( 1 );
    
    %% section 1. compute potential, including elestatic, and exchange potential
    
    
    % number of energy states below fermi level
    N = sum( epsilon_F >= E_schrodinger );
    % find normalization and |psi|^2*sqrt(m) term integrated over (x,y), normalized.
    [normalization, psi_sqrt_m] = normalize_and_sqrt_m_triangular( psi_schrodinger );
    % calculate the psi_sqrt_eps term needed to calculate the total
    % electron density
    psi_sqrt_eps_sum = psi_sqrt_eps_summation( N, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    % a row vector that calculates/contains the spatial regions of the
    % potential that are above the Fermi level; these spatial regions have
    % been ionized and contribute a positively-valued dopant density to the
    % Poisson equation
    heaviside_n_D = zeros( size(V_elestat') );
    heaviside_n_D( V_elestat > epsilon_F ) = 1;
    
    % SCALAR PDE EQUATION
    % f coefficient for the scalar PDE equation, solve the electrostatic
    % potentials, due to electrons, dopants, and polarizations. V_cb do not
    % need to be solved (assumed to remain unchanged through the iterations).
    f_poiss_1 = sprintf('n_D_func(x, y, %s) + n_e_func(x, y, %s)', mat2str(heaviside_n_D), mat2str(psi_sqrt_eps_sum));
    if debug
       disp("debug: iterate_run, showing the coefficients for the scalar poisson equation...");
       c_poiss
       a_poiss
       f_poiss_1
    end
    V_poiss = assempde(b_poiss, p_poiss, e_poiss, t_poiss, c_poiss, a_poiss, f_poiss_1);
    V_poiss = V_poiss / C;
    
    % N*N PDE EQUATIONS, for exchange potentials
    psi_sqrt_eps_xc_m = psi_sqrt_eps_xc( N, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    f_poiss_ex = f_ex_coefficient(N);
    b_poiss_ex = generate_dirichlet_boundary( N*N );
    if debug
       disp("debug: iterate_run, showing the coefficients for the N*N poisson equations for exchange potentials...");
       c_poiss
       a_poiss
       f_poiss_ex
    end
    % compute the exchange potentials
    V_xc = assempde(b_poiss_ex, p_poiss, e_poiss, t_poiss, c_poiss, a_poiss, f_poiss_ex);
    V_xc = V_xc / C;
    % update electrostatic potential, excluding exchange potentials
    V_elestat_old = V_elestat;
    V_elestat = V_poiss' + V_conduction_band(p_poiss(1, :), p_poiss(2, :));
    V_elestat = V_elestat - min(V_elestat);
    average_change = mean( abs(V_elestat-V_elestat_old) );
    
    %% section 2. solve schrodinger equations given the updated potentials
    % coefficient preparations for schordinger equations
    c_schrod = c_coefficient();
    % d coefficient is constant for all electrons and mesh points, a scalar
    d_schrod = 1;
    % build the a coefficient matrix
    a_schrod = a_coefficient_schrod(N, V_poiss, V_xc);
    r_schrod = [0, 1];
    b_schrod = generate_dirichlet_boundary( N );
    
    psi_schrodinger = []; % initialize the matrix psi_schrodinger
    E_schrodinger = []; % initialize the variable E_schrodinger
    epsilon_F = 0;
    v = [];
    l = [];
    
    while isempty(E_schrodinger) || E_schrodinger(end) <= epsilon_F
        % solve for N coupled schrodinger equations
        [v_add, l_add] = ...
                pdeeig(b_schrod, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod);
        
        % add newly-found electrons
        if ~isempty(l_add)
            
            v = [v, v_add];
            l = [l; l_add];
            
            if length(l) >= N
                if isempty(E_schrodinger)
                    % first, get N solutions and check if we've fulfilled
                    % charge neutrality constraint
                    [psi_schrodinger, E_schrodinger] = extract_solution(N, v, l);
                    [normalization, psi_sqrt_m] = normalize_and_sqrt_m_triangular(psi_schrodinger);
                    epsilon_F = find_epsilon_F( V_elestat, E_schrodinger, psi_sqrt_m );
                    % if not, we include more electrons in the solution if
                    % any
                    if E_schrodinger(end) > epsilon_F && length(l) > N
                        [psi_add, E_add] = update_solution(N, v_add(:, N+1:end), l(N+1:end));
                        E_schrodinger = [E_schrodinger; E_add];
                        psi_schrodinger = [psi_schrodinger, psi_add];
                    end
                else
                    % add new solutions
                    [psi_add, E_add] = update_solution(N, v_add, l_add);
                    E_schrodinger = [E_schrodinger; E_add];
                    psi_schrodinger = [psi_schrodinger, psi_add];
                end
                % calculate the normalization factor and the term of
                % |psi|^2*sqrt(m) integrated over (x,y)
                [normalization, psi_sqrt_m] = normalize_and_sqrt_m_triangular(psi_schrodinger);
                % calculate the fermi energy
                epsilon_F = find_epsilon_F( V_elestat, E_schrodinger, psi_sqrt_m );
                
            end
            
        end
        % increase the search range
        r_schrod(1) = r_schrod(2);
        r_schrod(2) = r_schrod(2) + 1;

    end
        
        
        
end