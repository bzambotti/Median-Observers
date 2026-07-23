function modes = generate_modes(q, p)
    % GENERATE_MODES Generates all possible delay configurations for the switched median observer.
    %
    % Inputs:
    %   q - Burst order (outliers can be up to q in a window of 2q+1)
    %   p - Number of output channels
    %
    % Outputs:
    %   modes - A matrix of size ( (2q+1)^p x p ) where each row is a valid delay configuration.
    %           Delays range from 0 to 2q.
    
    delays = 0:(2*q);
    
    if p == 1
        modes = delays(:);
        return;
    end
    
    % Create a cell array to hold the grid outputs
    grid_out = cell(1, p);
    
    % Generate the grid
    % We need to pass the same delay vector p times
    args = repmat({delays}, 1, p);
    [grid_out{:}] = ndgrid(args{:});
    
    % Flatten and concatenate to form the modes matrix
    modes = zeros(numel(grid_out{1}), p);
    for i = 1:p
        modes(:, i) = grid_out{i}(:);
    end
end
