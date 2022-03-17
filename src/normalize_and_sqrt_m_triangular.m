% a function to calculate the normalization pre-factor for all the
% wavefunctions and the |psi|^2*sqrt(m) term integrated over (x,y) used by
% find_epsilon_F.m to calculate the Fermi level.
function [normalization, psi_sqrt_m] = normalize_and_sqrt_m_triangular(psi_schrodinger)

    global debug;
    global t_schrod;
    global vector_of_masses;
    global triangle_indices_of_each_region_schrod;
    global triangle_areas_schrod;
    
    %%
    psi_squared = psi_schrodinger.^2;
    
    % a row vector of sqrt(mass) for each triangle
    triangle_sqrt_mass = zeros( size(triangle_areas_schrod) );
    for i = 1 : length(vector_of_masses)
        triangle_sqrt_mass( triangle_indices_of_each_region_schrod{i} ) = sqrt( vector_of_masses(i) );
    end

    normalization = zeros( size(psi_squared(1, :)) );

    psi_sqrt_m = zeros( size(psi_squared(1, :)) );

    for i = 1 : length(normalization)
        
        % a row vector of average |psi|^2 fir each triangle
        avg_psi_squared = 1/3 * ...
            ( psi_squared(t_schrod(1, :), i) + psi_squared(t_schrod(2, :), i) + psi_squared(t_schrod(3, :), i) );

        % matrix multiplication, the normalization factor for the i-th
        % wavefunction
        normalization(i) = 1 / sqrt( triangle_areas_schrod * avg_psi_squared );
        
        % |psi|^2*sqrt(m) term integrated over (x,y), normalized.
        psi_sqrt_m(i) = normalization(i)^2 * sum( triangle_areas_schrod.* avg_psi_squared'.* triangle_sqrt_mass );

    end
    
    %% debug purpose, showing the normalization and psi_sqrt_m
    if debug
        disp("debug: normalize_and_sqrt_m_triangular, showing the matrix of normalization and psi_sqrt_m below: ");
        normalization
        psi_sqrt_m
    end
    
    
end