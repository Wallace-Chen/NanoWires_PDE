% pick wavefunctions & energies from solutions of system of eigen-equations
function [psi, E] = extract_solution(N, v, l)
    
    % the number of electron states
    n_ele = min(N, length(l));
    % the number of mesh nodes
    Np = length( v(:,1) ) / N;
    % build mask matrix to speed up the solution extraction
    cells = cellmat(1, n_ele, Np, 1, 1);
    mask = blkdiag( cells{:} );
    
    % extract psi from eigenvectors v, pick the diagonal Np*1 block
    % matrices of v as wavefunctions:
    % psi_11 psi_12 psi_13
    % psi_21 psi_22 psi_23
    % psi_31 psi_32 psi_33
    % then chosen wavefunctions are [psi_11; psi_22; psi_33] with
    % eigen-energies increasing in order.
    psi = sum( v(1:Np*n_ele, 1:n_ele) .* mask, 2 );
    psi = reshape(psi, Np, n_ele);
    % slice the corresponding eigen-energies
    E = l(1:n_ele);
    
end