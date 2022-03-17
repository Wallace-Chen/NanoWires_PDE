% calculate the psi_sqrt_eps_ex term needed to calculate the source term of
% exchange potential:psi_i .* psi_j' * sqrt(epsilon_F - E_schrodinger_j) ).
% where psi_j' is the complex conjugate of the psi_j.
% @output, the output is a np * (N*N) matrix.
function mt = psi_sqrt_eps_xc( n_quanta, normalization, psi_schrodinger, epsilon_F, E_schrodinger )

    % use matrix operation to speed up
    % normalize the psi wavefunctions first
    psi_schrodinger = psi_schrodinger .* normalization;
    % slice the matrix
    psi_schrodinger = psi_schrodinger(:, 1:n_quanta);
    E_schrodinger = E_schrodinger(1 : n_quanta);
    
    mt = zeros( length(psi_schrodinger(:,1)), n_quanta*n_quanta );

    % exact exchange potential matrix, but not symmetric:
    %for i = 1 : n_quanta
    %   mt(:, (i-1)*n_quanta+1:i*n_quanta) = psi_schrodinger .* psi_schrodinger(:, i) * sqrt( epsilon_F - E_schrodinger(i) );
    %end
    
    % to make exchange potential matrix symmetric, we require V_{i,j} =
    % V_{j,i} = min(V_{i,j}, V_{j,i});
    for i = 1 : n_quanta
       mt(:, (i-1)*n_quanta+1:i*n_quanta) = psi_schrodinger .* psi_schrodinger(:, i) .* sqrt( epsilon_F - max( E_schrodinger(i), E_schrodinger' ) );
    end
    
end