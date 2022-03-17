% function to plot checkpoints variables
function plot_checkpoints(postfix)

    % global variables
    global ckp_n;
    global ckp_average_change;
    global ckp_epsilon_F;
    global ckp_E_schrodinger;
    
    global C;

    marker_size = 20;
    line_width = 3;
    label_size = 15;
    tick_size = 15;
    
    plotsFormat = {'MarkerSize', marker_size, 'LineWidth', line_width};
    labelFormat = {'fontsize', label_size};
    
    n_iter = length(ckp_n);
    x = 1 : n_iter;
    
    % plot number of occupied electrons
    figure
    plot(x, ckp_n, '.-', plotsFormat{:});
    xlabel('x (# of iter)', labelFormat{:});
    ylabel('y (# of electrons)', labelFormat{:});
    title('The number of occupied energies for iterations', labelFormat{:});
    set(gca,'fontsize', tick_size);
    grid on;
    saveas(gcf, ['./data/plots/n', postfix, '.jpeg']);
    
    % plot the average change of potentials
    figure
    plot(x(2:end), ckp_average_change(2:end) * C, '.-', plotsFormat{:});
    xlabel('x (# of iter)', labelFormat{:});
    ylabel('Average potential change (eV)', labelFormat{:});
    title('The average change of potentials (exchange potential excl.)', labelFormat{:});
    set(gca,'fontsize', tick_size);
    grid on;
    saveas(gcf, ['data/plots/ave_cha', postfix, '.jpeg'])
    
    % plot the fermi levels
    figure
    plot(x, ckp_epsilon_F, '.-', plotsFormat{:});
    xlabel('x (# of iter)', labelFormat{:});
    ylabel('Fermi Levels (1/C)', labelFormat{:});
    title('The fermi levels for iterations', labelFormat{:});
    set(gca,'fontsize', tick_size);
    grid on;
    saveas(gcf, ['data/plots/fer_lev', postfix, '.jpeg'])
    
    % plot occupied energies
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];
    ncolor = length(colors);
    figure
    for i = 1 : n_iter
        y = ckp_E_schrodinger{i};
        s = scatter(ones(1, length(y))*i, y, colors(mod(i,ncolor)+1), 'filled' ); 
        hold on;
    end
    xlabel('x (# of iter)', labelFormat{:});
    ylabel('Energy Levels (1/C)', labelFormat{:});
    title('Occupied energies for each iterations', labelFormat{:});
    set(gca,'fontsize', tick_size);
    grid on;
    saveas(gcf, ['data/plots/occ_ene', postfix, '.jpeg'])
    
end