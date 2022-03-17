% function to initiate the system, solve the indepedent schrodinger
% equation, find the proper number of electrons satisfying the charge
% neutrality equation.
function [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, V_xc, N] = initiate_system()

    global debug;
    
    global p_schrod; global e_schrod; global t_schrod;
    global p_poiss; global e_poiss; global t_poiss;
    global vector_of_eps;
    global number_of_interfaces;
    global vector_of_side_lengths;
    global C;
    global psi_sqrt_eps_xc_m;
    
    global save_file;
    
    disp("Initializing the system...");
    %% setup coefficient matrix
    % boundary matrix
    bmatrix = generate_dirichlet_boundary( 1 );
    % string representing c coefficient matrix
    % c_schrod = c_coefficient();
    c_schrod = 'c_coefficient_schrod(sd)';
    % string representing a coefficient matrix, only include conduction
    % band potential, no polarization, dopant, static or exchange
    % potentials.
    a_schrod = 'V_conduction_band(x,y)';
    % d coefficient is constant for all electrons and mesh points, a scalar
    d_schrod = 1;
    % initial guess for the potential, evaluated at mesh points
    V_elestat = V_conduction_band(p_poiss(1, :), p_poiss(2, :))';
    
    %% variable declaration
    psi_schrodinger = []; % initialize the matrix psi_schrodinger
    E_schrodinger = []; % initialize the variable E_schrodinger
    epsilon_F = 0;
    r_schrod = [0, 1];
    
    %% 1. solve independent schrodinger equation, keep adding electrons until the charge neutrality constraint satisfied.
    while isempty(E_schrodinger) || E_schrodinger(end) <= epsilon_F
        
        [~, psi_schrodinger_add, E_schrodinger_add] = ...
            evalc('pdeeig(bmatrix, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod)');
        
        % find new electrons with higher energies
        if ~isempty(E_schrodinger_add)
            
            psi_schrodinger = [psi_schrodinger, psi_schrodinger_add];
            E_schrodinger =   [E_schrodinger; E_schrodinger_add];
            
            % calculate the normalization factor and the term of
            % |psi|^2*sqrt(m) integrated over (x,y)
            [normalization, psi_sqrt_m] = normalize_and_sqrt_m_triangular(psi_schrodinger);
            % calculate the fermi energy
            epsilon_F = find_epsilon_F( V_elestat, E_schrodinger, psi_sqrt_m );
        end
        
        r_schrod(1) = r_schrod(2);
        r_schrod(2) = r_schrod(2) + 1;
        
    end
    
    % number of energy states below fermi level
    N = sum( epsilon_F >= E_schrodinger );
    
    % a row vector that calculates/contains the spatial regions of the
    % potential that are above the Fermi level; these spatial regions have
    % been ionized and contribute a positively-valued dopant density to the
    % Poisson equation
    heaviside_n_D = zeros( size(V_elestat') );
    heaviside_n_D( V_elestat > epsilon_F ) = 1;
    
    %% debug purpose, showing the wave functions, energies and fermi level
    if debug
       disp("Debug: initiate_system, the calculated eigen vectors, eigen-energies and fermi level are: ");
       psi_schrodinger
       E_schrodinger
       epsilon_F
    end
    
    %% 2. compute potential, including elestatic, and exchange potential
    
    [V_poiss, V_xc] = solve_potential(N, epsilon_F, E_schrodinger, psi_schrodinger, heaviside_n_D, 0.01);

    save(save_file, '-append');
    disp("Initialization finished.");
end