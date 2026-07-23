function d = generate_bursts(p, sim_steps, q, max_amp)
    % GENERATE_BURSTS Generates an outlier burst sequence of order q.
    % 
    % Inputs:
    %   p - Number of output channels
    %   sim_steps - Total simulation steps
    %   q - Outlier burst order (max q outliers in any 2q+1 window)
    %   max_amp - Maximum amplitude of the outliers
    %
    % Outputs:
    %   d - Disturbance matrix (p x sim_steps)

    d = zeros(p, sim_steps);
    N = 2 * q + 1;

    for j = 1:p
        k = 1;
        while k <= sim_steps
            % Decide whether to place a burst in this window
           % if rand() < 0.3 % 30% chance of a burst window
                % Number of outliers in this window (between 1 and q)
                num_outliers = randi([1, q]);

                % To strictly satisfy Definition 1, we must ensure that any sliding window of size N
                % has at most q outliers. We restrict the outliers to the first q indices of the block.
                % This ensures there are at least q+1 clear steps between bursts from adjacent blocks.
                end_idx = min(k + q - 1, sim_steps);
                window_len = end_idx - k + 1;

                if window_len > 0
                    actual_outliers = min(num_outliers, window_len);
                    positions = randperm(window_len, actual_outliers) + k - 1;

                    for pos = positions
                        % Random amplitude between -max_amp and max_amp
                        amp = (rand() * 2 - 1) * max_amp;
                        d(j, pos) = amp;
                    end
                end
          %  end
            % Move to the next window
            k = k + N;
        end
    end
end



