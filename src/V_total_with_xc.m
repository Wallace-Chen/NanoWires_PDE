% calculate the total potential including the exchange-correlation portion
%
% we need to make a separate routine called V_total_with_xc because we need
% to set the zero of the potential at the minimum of (V_poiss+V_xc) AFTER
% we sum them (setting the zero at the minimum of each separately and then
% adding them does not necessarily yield a total potential with its zero at
% the minimum)
function V_total_with_xc_val = V_total_with_xc(x, y, V_poiss)

global p_poiss
global t_poiss

% V_total_with_xc_val = interpolant_NaN(p_poiss, t_poiss, V_poiss, x, y)' +...
V_total_with_xc_val = griddata_NaN(p_poiss(1,:),p_poiss(2,:),V_poiss',x,y) +...
    V_conduction_band(x, y);

% sets the zero of the potential at the minimum of V_total_with_xc_val
V_total_with_xc_val = V_total_with_xc_val - min(V_total_with_xc_val, [], 'all');
