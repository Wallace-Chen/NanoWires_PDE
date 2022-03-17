% a function to build the a coefficient matrix for N coupled
% schrodinger equations with non-local exchange potentials included. 
% @output, a (N*N) * nr matrix. nr is the number of meshed triangles.
function a = a_coefficient_schrod(N, V_elestat, V_xc)

    global debug;
    
    global p_schrod; global t_schrod;
    global p_poiss; global t_poiss;
    
    global save_file;
    
    %%
    % get positions of schrod mesh triangles. a row vector
    p_mid_schrod = pdeintrp(p_schrod, t_schrod, p_schrod');
    x = p_mid_schrod(1, :);
    y = p_mid_schrod(2, :);
    
    % evaluate the conduction band potential at centers of triangles
    V_cb = V_conduction_band(x, y);
    % evaluate the elecstatic potentials at centers of triagnles, including
    % contributions from electrons and dopants.
    % V_elestat = interpolant_NaN(p_poiss, t_poiss, V_elestat, x, y)';
    V_elestat = griddata_NaN(p_poiss(1,:), p_poiss(2,:), V_elestat', x, y);
    
    % evaluate the exchange potentials at centers of triangles
    % V_xc_inter = interpolant_NaN(p_poiss, t_poiss, V_xc, x, y);
    % a = reshape(V_xc_inter, [], N*N)';
    np = length(p_poiss(1,:));
    a = ones(N*N, length(x));
    for i=1:N*N
       a( i, : ) = griddata_NaN(p_poiss(1,:), p_poiss(2,:), V_xc((i-1)*np+1 : i*np)', x, y);
    end
    
    % sum up all potentials
    n_col = 1;
    for i = 1 : N
       a(n_col, :) =  a(n_col, :) + V_cb + V_elestat;
       n_col = n_col + N + 1;
    end
    
    % set minimum potential point
    % important if we want to compare with the non-exchange potential case
    % we find the minimum of the diagonal terms of exchange potential
    % matrix, and set it zero
    % diag_indices = find(diag( ones(1, N) ));
    % mini = min( a(diag_indices, :), [], 'all' );
    
    % if the exchange potential matrix is made symmetric, we then can use a
    % N*(N+1)/2 column vector to represent the a coefficient matrix
    indices =  find( triu(reshape( 1:N*N, N, N )) );
    a = a( indices, : );
    
    % set zero-potential point
    a = a - min(a, [], 'all');
    % a = a - mini;
    
    save(save_file, '-append');
    
end