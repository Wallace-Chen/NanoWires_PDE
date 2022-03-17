% this function symmetrizes the coordinates using the symmetry of
% the core-shell nanowire
function [x, y] = symmetrize_coordinates(x, y)

    global debug;
    global degree_of_polygon;
    
    %%
    % place holder for the circle geometry
    if degree_of_polygon == 0
        x = x;
        y = y;
    % for the triangular geometry, only two-fold symmetry due to the
    % polarization terms
    elseif abs( degree_of_polygon ) == 3
        x = abs( x );
    elseif degree_of_polygon == 6
        % map to thee fourth quadrature
        x = abs( x );
        y = -abs( y );

        % reflection about a line:
        % https://math.stackexchange.com/questions/525082/reflection-across-a-line
        m= -1 / sqrt(3);
        reflection_matrix = 1/(m^2 + 1) * [1-m^2 2*m;2*m m^2-1]; 

        % find indices of coordinates that are above the y=-1/sqrt(3)*x line
        indices = find( y > m*x );
        % reflect about the y=-1/sqrt(3)*x line
        reflected_coordinates = reflection_matrix * [x(indices); y(indices)];

        x(indices) = reflected_coordinates(1, :);
        y(indices) = reflected_coordinates(2, :);

        m = -sqrt(3);
        reflection_matrix = 1/(m^2 + 1) * [1-m^2 2*m;2*m m^2-1];
        
        % find indices of coordinates that are above the y=-sqrt(3)*x line
        indices = find(y > m*x);

        % reflect about the y=-sqrt(3)*x line
        reflected_coordinates = reflection_matrix * [x(indices); y(indices)];

        x(indices) = reflected_coordinates(1, :);
        y(indices) = reflected_coordinates(2, :);
    end

end