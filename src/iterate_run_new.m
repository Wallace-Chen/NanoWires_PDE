% function to solve poisson and schrodinger equations iteratively until the
% potential gets converged.
% @param V_elestat: initial guess, contains only eletrcostatic potential.
function [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, V_xc, N] = ...
    iterate_run_new(N, V_poiss, V_xc)

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
    global damping_factor;
    
    global save_file;
    
    % make it global to have neat code
    global psi_sqrt_eps_xc_m;
    %%
    % build local potential
    V_elestat = V_poiss' + V_conduction_band(p_poiss(1, :), p_poiss(2, :));
    V_elestat = V_elestat - min(V_elestat);
    V_elestat = V_elestat';
    
    %% section 1. solve coupled schrodinger equations given the updated potentials
    fprintf('  1. Solving coupled schrodinger equations...\n');
    
    % coefficient preparations for schordinger equations
    % c_schrod = c_coefficient();
    c_schrod = 'c_coefficient_schrod(sd)';
    % d coefficient is constant for all electrons and mesh points, a scalar
    d_schrod = 1;
    n_iter = 1;
    
    input_N = N;
    input_V_poiss = V_poiss;
    n_increase = 1;
    
    scalar_E = [];
    scalar_psi = [];
    r_schrod_lower = 0;
    psi_schrodinger = [];
    E_schrodinger = [];
    epsilon_F = 0;
    while isempty(E_schrodinger) || E_schrodinger(end) <= epsilon_F
        %% 1.1 solve scalar eigen equation to find estimation of wavefunctions
        if n_iter > 1
            fprintf('    sub iter # %d, fermi level: %f, highest electron energy: %f\n', n_iter, epsilon_F, E_schrodinger(end));
            % now we solve scalar schrodinger poissons equations to find
            % these N electrons
            if isempty(scalar_E)
                [scalar_psi_add, scalar_E_add, r_schrod_lower] = solve_eig(n_increase + input_N, input_V_poiss, r_schrod_lower);
            else
                [scalar_psi_add, scalar_E_add, r_schrod_lower] = solve_eig(n_increase, input_V_poiss, r_schrod_lower);
            end
            scalar_psi = [scalar_psi, scalar_psi_add];
            scalar_E = [scalar_E; scalar_E_add];
            % below may be a poor treatment, since we don't update fermi
            % level here after we incorporate more electrons
            heaviside_n_D = zeros( size(V_elestat') );
            heaviside_n_D( V_elestat > epsilon_F ) = 1;
            N_old = N;
            N = length(scalar_E);
            V_xc_old = V_xc;
            V_poiss_old = V_poiss;
            [V_poiss, V_xc] = solve_potential(N, epsilon_F, scalar_E, scalar_psi, heaviside_n_D, 1);
            % update potential
            V_poiss = V_poiss_old + damping_factor * (V_poiss - V_poiss_old);
            V_xc = update_exchange_potential(V_xc_old, N_old, V_xc, N);
        end
        
        %% 1.2 solve multi-dimenstional eigen-equations
        fprintf('    Solving the multi-dimensional eigen-equations, input number of electrons: %d\n', N);
        % build the a coefficient matrix
        a_schrod = a_coefficient_schrod(N, V_poiss, V_xc);
        % a_schrod = 'V_total_with_xc(x,y)';
        b_schrod = generate_dirichlet_boundary( N );

        mini = 0;
        r_schrod = [mini, mini+1];

        psi_schrodinger = []; % initialize the matrix psi_schrodinger
        E_schrodinger = []; % initialize the variable E_schrodinger
        epsilon_F = 0;
        v = [];
        l = [];

        while isempty(E_schrodinger)
            %fprintf('      searching for solutiosn between [%d, %d], the current number of found electrons is: %d\n', r_schrod(1), r_schrod(2), length(l));
            % solve for N coupled schrodinger equations
            [~, v_add, l_add] = ...
                    evalc('pdeeig(b_schrod, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod)');

            % add newly-found electrons
            if ~isempty(l_add)

                v = [v, v_add];
                l = [l; l_add];

                if length(l) >= N
                    % get N solutions and check if we've fulfilled
                    % charge neutrality constraint
                    [psi_schrodinger, E_schrodinger] = extract_solution(N, v, l);
                    [normalization, psi_sqrt_m] = normalize_and_sqrt_m_triangular(psi_schrodinger);

                    epsilon_F = find_epsilon_F( V_elestat, E_schrodinger, psi_sqrt_m );
                end
            end
            % increase the search range
            r_schrod(1) = r_schrod(2);
            r_schrod(2) = r_schrod(2) + 1;
        end
        
        n_iter = n_iter + 1;
    end
    
        
    % number of energy states below fermi level
    N = sum( epsilon_F >= E_schrodinger );
    
    % a row vector that calculates/contains the spatial regions of the
    % potential that are above the Fermi level; these spatial regions have
    % been ionized and contribute a positively-valued dopant density to the
    % Poisson equation
    heaviside_n_D = zeros( size(V_elestat') );
    heaviside_n_D( V_elestat > epsilon_F ) = 1;
    
    
    %% section 2. compute potential, including elestatic, and exchange potential
    fprintf('  2. Solving poisson equations...\n');
    
    [V_poiss, V_xc] = solve_potential(N, epsilon_F, E_schrodinger, psi_schrodinger, heaviside_n_D, 1);
    
    save(save_file, '-append');
end