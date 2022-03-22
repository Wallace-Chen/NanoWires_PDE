
function [v, l] = debug_pdeeig(b_schrod, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod)
    interval = 5;

    v = [];
    l = [];
    r_schrod = [0, interval];

    while length(l) < 20
        fprintf('Searching between interval [%d, %d]\n', r_schrod(1), r_schrod(2));
        [~, v_add, l_add] = ...
                evalc('pdeeig(b_schrod, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod)');
            
        if ~isempty(l_add)
            v = [v, v_add];
            l = [l; l_add];
        end
        r_schrod(1) = r_schrod(2);
        r_schrod(2) = r_schrod(2) + interval;
    end

end