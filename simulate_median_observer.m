function [x_hat, e_norm] = simulate_median_observer(A, B, C, L, x, u, y, q)
    % SIMULATE_MEDIAN_OBSERVER Simulates the robust median-filtered observer.
    % 
    % Inputs:
    %   A, B, C - System matrices
    %   L - Observer gain(s). Can be a single matrix (fixed gain) or a cell array (mode-dependent gains).
    %   x - True state trajectory (n x sim_steps+1)
    %   u - Input trajectory (m x sim_steps)
    %   y - Measured output trajectory (p x sim_steps)
    %   q - Outlier burst order
    %
    % Outputs:
    %   x_hat - Estimated state trajectory (n x sim_steps+1)
    %   e_norm - Estimation error norm (1 x sim_steps+1)
    
    n = size(A, 1);
    p = size(C, 1);
    sim_steps = size(y, 2);
    N = 2 * q + 1;
    
    modes = generate_modes(q, p);
    num_modes = size(modes, 1);
    
    % Handle single gain vs mode-dependent gains
    if ~iscell(L)
        L_cell = cell(1, num_modes);
        for m = 1:num_modes
            L_cell{m} = L;
        end
    else
        L_cell = L;
    end
    
    x_hat = zeros(n, sim_steps + 1);
    e_norm = zeros(1, sim_steps + 1);
    
    % Initialize observer with some error
    x_hat(:, 1) = x(:, 1) + 2 * randn(n, 1);
    e_norm(1) = norm(x(:, 1) - x_hat(:, 1));
    
    % Memory buffer for the innovations (p x N)
    % Each row stores the past N innovations for an output channel
    % Index 1 is the most recent (k), index N is the oldest (k-2q)
    sigma = zeros(p, N);
    
    for k = 1:sim_steps
        % Compute current raw innovation
        tilde_y_k = y(:, k) - C * x_hat(:, k);
        
        % Shift buffer and insert new innovation
        sigma = [tilde_y_k, sigma(:, 1:N-1)];
        
        % Compute median for each channel and determine active delays
        median_val = zeros(p, 1);
        active_delays = zeros(1, p);
        
        for j = 1:p
            sorted_vals = sort(sigma(j, :));
            % Median is the (q+1)-th element since N = 2q+1
            med = sorted_vals(q + 1);
            median_val(j) = med;
            
            % Find the delay index (0 to 2q) that corresponds to this median value
            % If multiple match, pick the first one
            idx = find(sigma(j, :) == med, 1);
            active_delays(j) = idx - 1;
        end
        
        % Find the corresponding mode index
        mode_idx = 1;
        for m = 1:num_modes
            if all(modes(m, :) == active_delays)
                mode_idx = m;
                break;
            end
        end
        
        % Select the appropriate gain
        L_k = L_cell{mode_idx};
        
        % State update
        x_hat(:, k+1) = A * x_hat(:, k) + B * u(:, k) + L_k * median_val;
        
        % Compute error norm
        e_norm(k+1) = norm(x(:, k+1) - x_hat(:, k+1));
    end
end
