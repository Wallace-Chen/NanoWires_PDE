% function to plot checkpoints variables
function plot_debug(N, num, psi_schrodinger, E_schrodinger)

    %% global variables
    global l0;
    global p_schrod;
    global t_schrod;
    global C;
    
    global degree_of_polygon;

    %% parameter check
    if degree_of_polygon ~= 3 && degree_of_polygon ~= 6 && degree_of_polygon ~= 0
       disp("Error: generateMesh, the geometry not supported!");
       return;
    end
    
    %% plot all occupied wave-functions
    np = length(psi_schrodinger(:, 1)) / N;
%    [normalization, ~] = normalize_and_sqrt_m_triangular( psi_schrodinger );
    
    for j = 1 : num
        for i = 1 : N
            figure('visible','off');
            pdesurf(l0 * p_schrod, t_schrod, psi_schrodinger((i-1)*np+1:i*np,j) );
            view(0, 90);
            colormap jet;
            set(gcf, 'Renderer', 'zbuffer');
            xlabel('x (nm)');
            ylabel('y (nm)');
            title(['E_{',num2str(i),'-', num2str(j),'} = ',num2str(C*E_schrodinger(j)),' eV']);
            axis equal;
            axis off;

            saveas(gcf, ['./data/plots/debug/wavefunction_', num2str(i),'-', num2str(j), '.jpg']);
        end
    end
    
%     for j = 1 : num
%         [normalization, ~] = normalize_and_sqrt_m_triangular( reshape(psi_schrodinger(:,j), np, []) );
%         for i = 1 : N
%             figure('visible','off');
%             pdesurf(l0 * p_schrod, t_schrod, (normalization(i)*psi_schrodinger((i-1)*np+1:i*np,j)).^2 );
%             view(0, 90);
%             colormap jet;
%             set(gcf, 'Renderer', 'zbuffer');
%             xlabel('x (nm)');
%             ylabel('y (nm)');
%             title(['E_{',num2str(i),'-', num2str(j),'} = ',num2str(C*E_schrodinger(j)),' eV']);
%             axis equal;
%             axis off;
% 
%             saveas(gcf, ['./data/plots/debug/wavefunction_', num2str(i),'-', num2str(j), '_squared.jpg']);
%         end
%     end
        
    
end