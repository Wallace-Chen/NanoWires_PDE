% https://www.mathworks.com/help/pde/ug/pdeinterpolant.html
% https://www.mathworks.com/help/pde/ug/pde.pdeinterpolant.evaluate.html
% @param u: the name string of values of solution on the mesh (p, t), which
%           is the solution of PDEs;
function val = interpolant_NaN(p, t, u, x, y)

    F = pdeInterpolant(p, t, u);
    % https://www.mathworks.com/help/pde/ug/pde.pdeinterpolant.evaluate.html
    val = evaluate(F, x, y);
    val(isnan(val)) = 0;
    
end