% function to solve poisson and schrodinger equations iteratively until the
% potential gets converged.
% @param V_elestat: initial guess, contains only eletrcostatic potential.
function [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, V_xc, N] = ...
    iterate_run(N, V_poiss, V_xc)

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
    % build the a coefficient matrix
    a_schrod = a_coefficient_schrod(N, V_poiss, V_xc);
    % a_schrod = 'V_total_with_xc(x,y)';
    %% debug
    %% debug end
    
    % find the minimum of the whole potential matrix
    % mini = min(a_schrod, [], 'all');
    % fprintf('The minimum potential is: %f\n', mini);
    mini = 0;
    
    r_schrod = [mini, mini+10];
    b_schrod = generate_dirichlet_boundary( N );
    
    psi_schrodinger = []; % initialize the matrix psi_schrodinger
    E_schrodinger = []; % initialize the variable E_schrodinger
    epsilon_F = 0;
    v = [];
    l = [];
    
    n_iter = 1;
    while isempty(E_schrodinger) || E_schrodinger(end) <= epsilon_F
        % printing info
        ele_energy_high = 0;
        if ~isempty(E_schrodinger)
            ele_energy_high = E_schrodinger(end);
        end
%        fprintf('    %d iters: %d electrons, fermi energy: %f, the highest ele energy: %f, the higher end of range: %d \n',...
%            n_iter, length(E_schrodinger), epsilon_F, ele_energy_high, r_schrod(end));
        
        % solve for N coupled schrodinger equations
        [~, v_add, l_add] = ...
                evalc('pdeeig(b_schrod, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod)');
        
        %% debug purpose
%            fprintf('    number of newly-found electrons is: %d\n', length(l_add));
%            fprintf('    number of found electrons is: %d\n', length(l));
        %%
            
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
                    
                    save(save_file, '-append');
                    
                    epsilon_F = find_epsilon_F( V_elestat, E_schrodinger, psi_sqrt_m );
                    % if not, we include more electrons in the solution if
                    % any
                    if E_schrodinger(end) < epsilon_F && length(l) > N
                        [psi_add, E_add] = update_solution(N, v(:, N+1:end), l(N+1:end));
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
        r_schrod(2) = r_schrod(2) + 10;
        
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