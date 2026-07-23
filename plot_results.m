function plot_results(x, x_hat, e_norm, y, sim_steps, active_delays_log)
    % PLOT_RESULTS Visualizes the simulation results.
    % 
    % Inputs:
    %   x - True state trajectory (n x sim_steps+1)
    %   x_hat - Estimated state trajectory (n x sim_steps+1)
    %   e_norm - Estimation error norm (1 x sim_steps+1)
    %   y - Measured output trajectory (p x sim_steps)
    %   sim_steps - Total simulation steps
    %   active_delays_log - Log of active delays per channel (p x sim_steps) (Optional)

    if nargin < 6
        active_delays_log = [];
    end
    
    n = size(x, 1);
    p = size(y, 1);
    time = 0:sim_steps;
    
    figure('Position', [100, 100, 1200, 1000]);

    num_subplots = 3;
    if ~isempty(active_delays_log)
        num_subplots = 4;
    end
    
    % Plot Error Norm (Log Scale)
    subplot(num_subplots, 1, 1);
    semilogy(time, e_norm, 'LineWidth', 2);
    title('Estimation Error Norm ||e_k||');
    xlabel('Time Step k');
    ylabel('Norm (log scale)');
    grid on;
    
    % Plot States
    subplot(num_subplots, 1, 2);
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
    subplot(num_subplots, 1, 3);
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
    
    % Plot Active Delays
    if ~isempty(active_delays_log)
        subplot(num_subplots, 1, 4);
        hold on;
        for j = 1:p
            stairs(1:sim_steps, active_delays_log(j, :), 'Color', colors(j, :), 'LineWidth', 1.5, 'DisplayName', sprintf('Channel %d', j));
        end
        title('Median Filter Active Delay Selection Index');
        xlabel('Time Step k');
        ylabel('Delay Index (0 to 2q)');
        yticks(0:(max(active_delays_log(:))+1));
        legend('Location', 'best');
        grid on;
        hold off;
    end

    % Adjust layout
    try
        sgtitle('Median-Filtered Observer Performance');
    catch
        % Octave might not support sgtitle
    end
end
