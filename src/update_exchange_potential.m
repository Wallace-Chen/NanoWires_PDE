% implement the update of the exchange potential matrix
function V_xc_updated = update_exchange_potential(V_xc_old, N_old, V_xc_new, N_new)
    
    global damping_factor;
    
    % reshape a np*N*N column vector into a N*N by np matrix for later use 
    V_xc_old = reshape(V_xc_old, [], N_old^2)';
    V_xc_new = reshape(V_xc_new, [], N_new^2)';
    
    V_xc = V_xc_old;
    % if the dimension get reduced, we trim the old matrix
    if N_old > N_new
        indices = reshape( 1:N_old^2, N_old, N_old );
        indices(N_new+1:end, :) = 0;
        indices(:, N_new+1:end) = 0;
        V_xc = V_xc_old( find(indices), : );
    end
    % we expand the old matrix if its dimension is less than that of new
    % matrix
    if N_old < N_new
        N_diff =  N_new - N_old;
        indices = reshape( 1:N_new^2, N_new, N_new );
        indices(N_old+1:end, :) = 0;
        indices(:, N_old+1:end) = 0;
        
        V_xc = zeros( size(V_xc_new, 1), size(V_xc_new, 2) );
        V_xc( find(indices), : ) = V_xc_old;
        indices(N_old+1:end, :) = -1;
        V_xc( find(indices==0) ,: ) = repmat(V_xc_old( end+1-N_old:end , : ),  N_diff, 1);
        indices(N_old, :) = -2;
        V_xc( find(indices==-1),: ) = repelem( V_xc(find(indices==-2), :), N_diff, 1);
    end
    % now the dimension of V_xc should be equal to V_xc_new
    V_xc_updated = V_xc + damping_factor * (V_xc_new - V_xc);
    
    % reshape back to column vector
    V_xc_updated = reshape( V_xc_updated', [], 1 );
    
end