function [L_m, P, N_slack, feasible] = codesign_switched(A, C, q)
    % CODESIGN_SWITCHED Synthesizes mode-dependent observer gains using the exact Switched framework.
    % 
    % Inputs:
    %   A - System matrix (n x n)
    %   C - Output matrix (p x n)
    %   q - Outlier burst order (scalar)
    %
    % Outputs:
    %   L_m - Cell array of synthesized observer gains (n x p) for each mode
    %   P - Common Lyapunov matrix (Nn x Nn)
    %   N_slack - Common slack matrix
    %   feasible - Boolean indicating if the LMIs were solved successfully
    
    n = size(A, 1);
    p = size(C, 1);
    N = 2 * q + 1;
    
    % Generate all possible modes (delay combinations)
    modes = generate_modes(q, p);
    num_modes = size(modes, 1);
    
    % Define \bar{A}, \bar{B}
    A_bar = zeros(N*n, N*n);
    A_bar(1:n, 1:n) = A;
    for i = 1:(N-1)
        A_bar(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = eye(n);
    end
    
    B_bar = zeros(N*n, n);
    B_bar(1:n, 1:n) = eye(n);
    
    % YALMIP variables
    P_var = sdpvar(N*n, N*n, 'symmetric');
    N1 = sdpvar(n, n, 'full');
    N12 = sdpvar(n, 2*q*n, 'full');
    N2 = sdpvar(2*q*n, 2*q*n, 'full');
    N_var = [N1, N12; zeros(2*q*n, n), N2];
    
    Y_vars = cell(1, num_modes);
    for m = 1:num_modes
        Y_vars{m} = sdpvar(n, p, 'full');
    end
    
    Constraints = [P_var >= 1e-6*eye(N*n)];
    
    % Build LMI for each mode
    for m = 1:num_modes
        delay_cfg = modes(m, :);
        
        % Construct \bar{C}_m for this mode
        C_bar_m = zeros(p, N*n);
        for j = 1:p
            delay = delay_cfg(j);
            % Place C(j, :) at the corresponding delay block (delay 0 is block 1, delay d is block d+1)
            C_bar_m(j, delay*n + 1 : (delay+1)*n) = C(j, :);
        end
        
        % LMI components
        M11 = -P_var;
        M12 = A_bar' * N_var' - C_bar_m' * Y_vars{m}' * B_bar';
        M22 = P_var - N_var - N_var';
        
        LMI_m = [M11, M12;
                 M12', M22];
                 
        Constraints = [Constraints, LMI_m <= -1e-6*eye(size(LMI_m))];
    end
    
    % Optimization
    options = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(Constraints, [], options);
    
    if sol.problem == 0
        feasible = true;
        P = value(P_var);
        N_slack = value(N_var);
        
        L_m = cell(1, num_modes);
        for m = 1:num_modes
            L_m{m} = value(N1) \ value(Y_vars{m});
        end
    else
        feasible = false;
        P = NaN(N*n, N*n);
        N_slack = NaN(N*n, N*n);
        L_m = cell(1, num_modes);
    end
end
