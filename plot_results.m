function plot_results(x, x_hat, e_norm, y, sim_steps)
    % PLOT_RESULTS Visualizes the simulation results.
    % 
    % Inputs:
    %   x - True state trajectory (n x sim_steps+1)
    %   x_hat - Estimated state trajectory (n x sim_steps+1)
    %   e_norm - Estimation error norm (1 x sim_steps+1)
    %   y - Measured output trajectory (p x sim_steps)
    %   sim_steps - Total simulation steps
    
    n = size(x, 1);
    p = size(y, 1);
    time = 0:sim_steps;
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot Error Norm (Log Scale)
    subplot(3, 1, 1);
    semilogy(time, e_norm, 'LineWidth', 2);
    title('Estimation Error Norm ||e_k||');
    xlabel('Time Step k');
    ylabel('Norm (log scale)');
    grid on;
    
    % Plot States
    subplot(3, 1, 2);
    hold on;
    colors = lines(n);
    for i = 1:n
        plot(time, x(i, :), '-', 'Color', colors(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('x_%d', i));
        plot(time, x_hat(i, :), '--', 'Color', colors(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('x\\_hat_%d', i));
    end
    title('True vs Estimated States');
    xlabel('Time Step k');
    ylabel('State Value');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Plot Outputs with Outliers
    subplot(3, 1, 3);
    hold on;
    colors = lines(p);
    for j = 1:p
        plot(1:sim_steps, y(j, :), 'x', 'Color', colors(j, :), 'MarkerSize', 6, 'DisplayName', sprintf('y_%d', j));
    end
    title('Measured Outputs (with outlier bursts)');
    xlabel('Time Step k');
    ylabel('Output Value');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Adjust layout
    try
        sgtitle('Median-Filtered Observer Performance');
    catch
        % Octave might not support sgtitle
    end
end
