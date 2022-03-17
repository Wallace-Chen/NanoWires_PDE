% a function to return a string representing a coefficient matrix for the
% PDE equation, a contains all potentials experienced by the electron,
% including v_cb, v_stat and v_ex. v_cb (conduction band) is constant for
% all eletrcons, v_stat (eletrostatic potential due to dopants, electrons and
% polatization) is different for each electron, can be represented by a
% diagonal matrix. v_ex is the coupling term from the exchange potential,
% represented by a symmetrical matrix.

% @param v_cb:   a string representing the v_cb, which is the same for all
%                equation. Empty strings means no such potential.
% @param v_stat: string of Np*N column vector, the first Np elements represent the
%                electrostatic potential by 
%                the first electron. Np: the number of the mesh points. N:
%                number of system. empty means no potential.
%                Note v_stat is the name (string) of your column vector,
%                not the matrix itself!! 
% @param v_xc:   string of Np*(N*(N+1)/2) column vector, the first Np elements
%                represent the exchange potential between two electrons
%                1,1. Np: the number of the mesh points. consider 3 eles 
%                and the matrix:
%                v11  v12  v13
%                v21  v22  v23
%                v31  v32  v33
%                then your input v_xc should be: [v11; v12; v22; v13; v23;
%                v33]. each vij is N-sized column vector. 
%                empty means no potential.
%                Note v_xc is the name (string) of your column vector, not
%                the matrix itself!!
% @param p:      matrix of mesh points to be used;
% @param t:      matrix of mesh triangles to be used.
% output: a string representing the total potential.
%         if total potential is constant (in terms of different electrons),
%         then a scalar is returned;
%         if the total potential has no coupling terms, then a N-sized
%         column vector is returned;
%         if the total potential has symmetric coupling terms (always true
%         in this program), a N*(N+1)/2 sized column vector is returned.
function a = a_coefficient(N, v_cb, v_stat, v_xc, p, t)

    global debug;
    %% parameter pre-processing
    
    % test if v_cb is empty char or string.
    cb_empty = isempty(v_cb) || v_cb == "";
    % test if v_stat is empty char, string.
    stat_empty = isempty(v_stat) || v_stat == "";
    % test if v_xc is empty char, string.
    xc_empty = isempty(v_xc) || v_xc == "";
    a = '';
    
    p = mat2str(p);
    t = mat2str(t);
    %% build the coefficient matrix
    % a (N*(N+1)/2) column vector
    if ~xc_empty
        a = build_matrix_coefficient(N*(N+1)/2, v_xc, p, t);
    end
    if ~stat_empty
       % a N-sized column vector
       tmp = build_matrix_coefficient(N, v_stat, p, t);
       % has to combine (N*(N+1)/2) column vector with N-sized column
       % vector.
       if ~isempty(a)
           idx = 0;
           for i = 1 : N
               idx = idx + i;
               a(idx) = a(idx) + " + " + tmp(i);
           end
       else
           a = tmp;
       end
    end
    if ~cb_empty
        if isempty(a)
            a = v_cb;
        else
            idx = 0;
            for i = 1 : N
                if length(a) == N
                    idx = i;
                else
                    idx = idx + i;
                end
                a(idx) = a(idx) + " + " + v_cb;
            end
        end
    end
    
    
    if isempty(a) || a == ""
       disp('Warning: a_coefficient, a coefficient is 0, make sure it is expected.');
       a = 0;
    else
        a = char( convertStringsToChars(a) );
    end
    
    %% debug purpose: showing a coefficient
    if debug
        disp("Debug: a_coefficient, a coefficient is: ...");
        a
    end
    

end
