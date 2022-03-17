% when length(l) > N, add eigenvector of last equation of each solutions 
function [psi_add, E_add] = update_solution(N, v, l)

    Np = length( v(:,1) ) / N;
    
    E_add = l;
    psi_add = v(end-Np+1:end, :);

end