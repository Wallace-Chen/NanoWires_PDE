% function will generate dirichlet boundary matrix for the geometry with
% degree_of_polygon vertices, the dimenstion of system is given by N
function bmatrix = generate_dirichlet_boundary( N )

    global debug;

    global degree_of_polygon;
    %% parameter processing
    % the number of dirichlet conditions
    M = N;
    nrow = 1;
    % the number of edgens for the polygon
    nedges = degree_of_polygon;
    if degree_of_polygon == 0
       nedges = 4;
    end
    
    %% build the boundary matrix
    
    bmatrix = ones(2+2*(N^2+N+N*M+M), nedges);
    % the dimension of the system
    bmatrix(nrow, :) = ones(1, nedges) * N;
    nrow = nrow + 1;
    % the number of dirichlet boundary conditions
    bmatrix(nrow, :) = ones(1, nedges) * M;
    nrow = nrow + 1;
    % length for the strings representing q
    bmatrix(nrow:nrow+N^2-1, :) = ones(N^2, nedges);
    nrow = nrow + N^2;
    % length for the strings representing g
    bmatrix(nrow:nrow+N-1, :) = ones(N, nedges);
    nrow = nrow + N;
    % length for the strings representing h
    bmatrix(nrow:nrow+M*N-1, :) = ones(N*M, nedges);
    nrow = nrow + N*M;
    % length for the strings representing r
    bmatrix(nrow:nrow+M-1, :) = ones(M, nedges);
    nrow = nrow + M;
    
    % the expression for q
    bmatrix(nrow:nrow+N^2-1, :) = ones(N^2, nedges) * 48;
    nrow = nrow + N^2;
    % the expression for g
    bmatrix(nrow:nrow+N-1, :) = ones(N, nedges) * 48;
    nrow = nrow + N;
    % the expression for h
    bmatrix(nrow:nrow+M*N-1, :) = repmat( reshape( diag(ones(1, N)) + 48, [], 1), 1, nedges);
    nrow = nrow + N*M;
    % length for the strings representing r
    bmatrix(nrow:nrow+M-1, :) = ones(M, nedges) * 48;
    
    
    %% debug: showing the boundary matrix
    if debug
        disp("Debug: generateDirichletBoundary, showing the boundary matrix...");
        bmatrix
    end
    
end