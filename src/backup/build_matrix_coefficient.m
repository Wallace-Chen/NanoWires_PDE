% a helper to return the coefficient matrix when it is diagonal.
% the function will automatically add V_conduction_band(x,y).

% @param vs, string of the matrix of the N*np column vector, which is the
%            solution of N PDEs. np is the number of mesh points;
% @param p: string of matrix of mesh points to be used;
% @param t: string of matrix of mesh triangles to be used
function s = build_matrix_coefficient(N, vs, p, t)

    s = strings(N, 1);
    
    for i = 1 : N
       s(i) = sprintf("interpolant_NaN(%s, %s, %s, %d, x, y)", p, t, vs, i);
    end
        
    % convert to N-sized char arrays
    % s = char( convertStringsToChars(s) );

end