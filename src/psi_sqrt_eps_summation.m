% calculate the psi_sqrt_eps term needed to calculate the total electron
% density: sum( |psi|^2 * sqrt(fpsilon_F - E_schrodinger) )
function summation = psi_sqrt_eps_summation( n_quanta, normalization, psi_schrodinger, epsilon_F, E_schrodinger )

    % normalize the psi wavefunctions first
    psi_schrodinger = psi_schrodinger .* normalization;

    % use matrix operation instead
    psi_schrodinger = psi_schrodinger(:, 1:n_quanta);
    E_schrodinger = E_schrodinger(1 : n_quanta);
    psi_schrodinger = psi_schrodinger .^ 2 .* sqrt( epsilon_F-E_schrodinger' );
    summation = sum(psi_schrodinger, 2); % sum along row axis
    
%     for i = 1 : n_quanta
% 
%         % the 1/pi*sqrt_mass prefactor is calculated later in the function
%         % n_e_prefactor.m
%         summation = summation + psi_schrodinger(:,i) .^2 * sqrt( epsilon_F-E_schrodinger(i) );
%     end

end