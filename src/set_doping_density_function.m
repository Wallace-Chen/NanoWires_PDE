% this function sets the doping density function to be used in the HADOKEN
% program. The default value given below will generate the plots in Figure
% 5 of the main text when main_scf_neumann is executed
function n_D = set_doping_density_function(x, y)

    % doping density (in units of 10^18 cm^-3, or 10^24/m^3). This has to be
    % inputted by the user as an array-smart expression. For example, if you
    % want to specify a doping density of the form 0.4*exp(-0.1*(x^2+y^2)), the
    % above line would be changed to n_D=0.4*exp(-0.1*(x.^2+y.^2));
    n_D = 0.2;


    if length( n_D ) == 1
        initialize_matrix = zeros( size(x) );
        initialize_matrix(:) = n_D;
        n_D = initialize_matrix;
    end
    % the commands above do not need to be modified by the user; they just
    % ensure that n_D has the correct dimensions if a scalar value is used for
    % the variable n_D (n_D must have the same size as the number of grid
    % points specified by the variable "x")
    
end