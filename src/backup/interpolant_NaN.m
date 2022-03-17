% https://www.mathworks.com/help/pde/ug/pdeinterpolant.html
% https://www.mathworks.com/help/pde/ug/pde.pdeinterpolant.evaluate.html
% @param u: the name string of values of solution on the mesh (p, t), which
%           is the solution of PDEs;
% @param idx: indicate the index of equation's solution to interpolate and evaluate.
function val = interpolant_NaN(p, t, u, idx, x, y)

    np = length(p(1, :));
    F = pdeInterpolant(p, t, u( (idx-1)*np+1 : idx*np ));
    val = evaluate(F, x, y)';
    val(isnan(val)) = 0;
    
end