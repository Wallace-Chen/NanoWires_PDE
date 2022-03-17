% function to create PDE model and load it with user-defined geometry
function generateMesh(degree_of_polygon, vector_of_side_lengths)

    % for debug purpose in development phase, set to false in production
    % phase 
    debug = true;

    %% parameter check
    if degree_of_polygon ~= 3 & degree_of_polygon ~= 6 & degree_of_polygon ~= 0
       disp("Error: generateMesh, the geometry not supported!");
       return;
    end
    
    %% initiate geometry matrix
    % get the number of layers
    number_of_interfaces = length(vector_of_side_lengths);
    
    % the number of rows of geometry description matrix
    nrows = 2 + 2*degree_of_polygon;
    % special treatment for the circle
    if degree_of_polygon == 0
        nrows = 4;
    end
    
    % the geometry decription matrix
    geom = zeros(nrows, number_of_interfaces);
    
    ns = strings([number_of_interfaces, 1]);
    sf = 'S1';
    %%  build the geometry matrix
    % construct the geometry description matrix:
    % https://www.mathworks.com/help/pde/ug/decsg.html 
    for i = 1:number_of_interfaces
       
        ns(i) = sprintf("S%d", i);
        if i > 1
            sf = [sf, ' + ', char(ns(i))];
        end
        edge_length = vector_of_side_lengths(i);
        % treatment for the circle
        if degree_of_polygon == 0
            geom(1, i) = 1;
            geom(2, i) = 0;
            geom(3, i) = 0;
            geom(4, i) = edge_length;
            continue;
        end
       
        % general approach for the polygon
        geom(1, i) = 2; geom(2, i) = degree_of_polygon;
        % treatment for triangle
        % TODO! assumes a N-face configuation for now, need add Ga-face
        if degree_of_polygon == 3 
            geom(3, i) = 0;             geom(6, i) = sqrt(3)/3 * edge_length;
            geom(4, i) = edge_length/2; geom(7, i) =-sqrt(3)/6 * edge_length;
            geom(5, i) =-edge_length/2; geom(8, i) =-sqrt(3)/6 * edge_length;
        % treatment for hexagon
        elseif degree_of_polygon == 6
            geom(3, i) =-edge_length/2; geom(9, i) = sqrt(3)/2 * edge_length;
            geom(4, i) = edge_length/2; geom(10,i) = sqrt(3)/2 * edge_length;
            geom(5, i) = edge_length;   geom(11,i) = 0;
            geom(6, i) = edge_length/2; geom(12,i) =-sqrt(3)/2 * edge_length;
            geom(7, i) =-edge_length/2; geom(13,i) =-sqrt(3)/2 * edge_length;
            geom(8, i) =-edge_length;   geom(14,i) = 0;
        end               
    end

    ns = char(ns)';
    
    % create geometry
    gd = decsg(geom, sf, ns);
    
    %% debug purpose
    if debug
        disp('Debug: generateMesh, showing the geometry...');
        gd
        pdegplot(gd, 'EdgeLabels', 'On', ...
                     'FaceLabels', 'On', ...
                     'VertexLabels', 'On');
    end
    
end