% function to create user-defined geometry and triangular meshes
function [p, e, t, triangle_indices_of_each_region, triangle_areas] = generate_mesh()

    global debug;
    
    global degree_of_polygon;
    global vector_of_side_lengths;
    global number_of_triangles;
    
    %% parameter check
    if degree_of_polygon ~= 3 && degree_of_polygon ~= 6 && degree_of_polygon ~= 0
       disp("Error: generateMesh, the geometry not supported!");
       return;
    end
    
    %% initiate decomposed geometry matrix
    % get the number of layers
    number_of_interfaces = length(vector_of_side_lengths);
    outer_area = 3 * sqrt(3) * vector_of_side_lengths(1) ^ 2 / 2;
    
    % the number of rows of geometry description matrix
    nrows = 7;
    nedges = degree_of_polygon;
    % special treatment for the circle
    if degree_of_polygon == 3
        outer_area = qrt(3) * vector_of_side_lengths(1) ^ 2 / 4;
    elseif degree_of_polygon == 0
        nrows = 10;
        nedges = 4;
        outer_area = pi * vector_of_side_lengths(1) ^ 2;
    end
    
    % the geometry decription matrix
    geom = zeros(nrows, number_of_interfaces * nedges);
    
    %%  build the decomposed geometry matrix
    % construct the geometry description matrix:
    % https://www.mathworks.com/help/pde/ug/decsg.html#bu_fft3-3
    for i = 1:number_of_interfaces
        j = (i-1) * nedges + 1;
        edge_length = vector_of_side_lengths(i);
        interface_coordinates_x = [];
        interface_coordinates_y = [];
        % treatment for triangle
        % TODO! assumes a N-face configuation for now, need add Ga-face
        if degree_of_polygon == 3 
            geom(1, j:j+nedges-1) = ones(1, nedges) * 2;
            interface_coordinates_x = [0 -edge_length/2 edge_length/2];
            interface_coordinates_y = [sqrt(3)/3 * edge_length -sqrt(3)/6 * edge_length -sqrt(3)/6 * edge_length];
        % treatment for hexagon
        elseif degree_of_polygon == 6
            geom(1, j:j+nedges-1) = ones(1, nedges) * 2;
            interface_coordinates_x = [-edge_length/2 -edge_length -edge_length/2 edge_length/2 edge_length edge_length/2];
            interface_coordinates_y = [sqrt(3)/2 * edge_length 0 -sqrt(3)/2 * edge_length -sqrt(3)/2 * edge_length 0 sqrt(3)/2 * edge_length];
        % treatment for circle
        elseif degree_of_polygon == 0
            interface_coordinates_x = [0 -edge_length 0 edge_length];
            interface_coordinates_y = [edge_length 0 -edge_length 0];
            geom(1, j:j+nedges-1) = ones(1, nedges);
            geom(8, j:j+nedges-1) = zeros(1, nedges);
            geom(9, j:j+nedges-1) = zeros(1, nedges);
            geom(10, j:j+nedges-1) = ones(1, nedges) * edge_length;
        end
        % starting x-coordinate
        geom(2, j:j+nedges-1) = interface_coordinates_x;
        % ending x-coordinate
        geom(3, j:j+nedges-1) = interface_coordinates_x([2:end 1]);
        % starting y-coordinate
        geom(4, j:j+nedges-1) = interface_coordinates_y;
        % ending y-coordinate
        geom(5, j:j+nedges-1) = interface_coordinates_y([2:end 1]);
        % left minimal region label
        geom(6, j:j+nedges-1) = ones(1, nedges) * i;
        % right minimal region label
        geom(7, j:j+nedges-1) = ones(1, nedges) * (i-1);
    end
    
    %% debug purpose
    g=zeros(7,number_of_interfaces*6);

    g(1,:)=2;
    column_index=1;

    for j=1:number_of_interfaces

        interface_side_length=vector_of_side_lengths(j);
        % extract out the side length of the outermost interface first and
        % "works inwards" towards the innermost core region

        interface_coordinates_x=[-interface_side_length/2 interface_side_length/2 interface_side_length interface_side_length/2 -interface_side_length/2 -interface_side_length];
        interface_coordinates_y=[sqrt(3)/2*interface_side_length sqrt(3)/2*interface_side_length 0 -sqrt(3)/2*interface_side_length -sqrt(3)/2*interface_side_length 0];

        g(2,column_index:column_index+5)=interface_coordinates_x;

        g(3,column_index:column_index+5)=interface_coordinates_x([2:end 1]);

        g(4,column_index:column_index+5)=interface_coordinates_y;

        g(5,column_index:column_index+5)=interface_coordinates_y([2:end 1]);

        g(6,column_index:column_index+5)=(j-1)*ones(1,6);

        g(7,column_index:column_index+5)=j*ones(1,6);

        column_index=column_index+6;

    end
    
    geom = g;
    %% debug end
    
    %% debug: showing geometry
    if debug
        disp('Debug: generateMesh, showing the geometry...');
        geom
        figure;
        pdegplot(geom, 'EdgeLabels', 'On', ...
                     'FaceLabels', 'On', ...
                     'VertexLabels', 'On');
    end
    %% debug end
    
    area_of_individual_triangle_mesh = outer_area / number_of_triangles;
    edge_size = sqrt(4 * area_of_individual_triangle_mesh / sqrt(3));
    
    % create 2D triangular mesh,
    % https://www.mathworks.com/help/pde/ug/initmesh.html 
    [p, e, t] = initmesh(geom, 'Hmax', edge_size, 'Jiggle', 'on');
    
    %% the following block code may not be necessary...
    edge_matrix_indices = ['e(5,:)==', num2str(nedges+1)];
    for j = nedges + 2 : number_of_interfaces * nedges
        edge_matrix_indices = [edge_matrix_indices,' | e(5,:)==',num2str(j)];
    end
    % the edge matrix includes geometric information for all the inner
    % interfaces so that initmesh will include finite element grid points on
    % these edges (i.e., this allows us to resolve the abrupt interfaces more
    % finely). However, we do not want to set boundary conditions on these
    % inner interfaces so we eliminate points in the edge matrix that
    % correspond to all the inner core boundaries
    e(:, eval(edge_matrix_indices))=[];
    e(7, :)=1;
    
    %% debug: showing the mesh
    if debug
        disp("Debug: generateMesh, edge_size is " + num2str(edge_size));
        disp("Debug: generateMesh, size of p is: "); size(p)
        disp("Debug: generateMesh, size of e is: "); size(e)
        disp("Debug: generateMesh, size of t is: "); size(t)
        disp('Debug: generateMesh, showing the triangle mesh...');
        figure;
        pdemesh(10 * p, e, t);
        xlabel('x (nm)');
        ylabel('y (nm)');
        axis equal;
    end
    %% debug end
 
    for i = 1:number_of_interfaces

        % find out indices of triangles that reside within various regions in
        % the multi-coreshell nanowire. This array is ordered as follows: the
        % first element contains indices of triangles within the outermost
        % region, and the subsequent elements contains indices of triangles as
        % one "works inwards" towards the innermost core region
        triangle_indices_of_each_region{i} = find( t(4, :) == i );

    end
    
    % t(4, :)=1;
    
    x_A = p(1, t(1, :));
    x_B = p(1, t(2, :));
    x_C = p(1, t(3, :));
    y_A = p(2, t(1, :));
    y_B = p(2, t(2, :));
    y_C = p(2, t(3, :));
    
    % a row vector of triangle masses
    triangle_areas = 1/2 * abs( (x_C-x_A).* (y_B-y_A) - (x_B-x_A).* (y_C-y_A) );
    
    
end