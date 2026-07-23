function [P, feasible] = check_switched_analysis(A, C, L, q)
    % CHECK_SWITCHED_ANALYSIS Verifies stability of given observer gain(s) under arbitrary switching.
    % 
    % Inputs:
    %   A - System matrix (n x n)
    %   C - Output matrix (p x n)
    %   L - Observer gain(s). Can be a single matrix (n x p) for fixed gain, 
    %       or a cell array of matrices for mode-dependent gains.
    %   q - Outlier burst order (scalar)
    %
    % Outputs:
    %   P - Common Lyapunov matrix (Nn x Nn)
    %   feasible - Boolean indicating if the LMI was solved successfully
    
    n = size(A, 1);
    p = size(C, 1);
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
    
    % Define \bar{A}, \bar{B}
    A_bar = zeros(N*n, N*n);
    A_bar(1:n, 1:n) = A;
    for i = 1:(N-1)
        A_bar(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = eye(n);
    end
    
    B_bar = zeros(N*n, n);
    B_bar(1:n, 1:n) = eye(n);
    
    % YALMIP variable
    P_var = sdpvar(N*n, N*n, 'symmetric');
    Constraints = [P_var >= 1e-6*eye(N*n)];
    
    % Build LMI for each mode
    for m = 1:num_modes
        delay_cfg = modes(m, :);
        
        C_bar_m = zeros(p, N*n);
        for j = 1:p
            delay = delay_cfg(j);
            C_bar_m(j, delay*n + 1 : (delay+1)*n) = C(j, :);
        end
        
        A_cl_m = A_bar - B_bar * L_cell{m} * C_bar_m;
        
        LMI_m = A_cl_m' * P_var * A_cl_m - P_var;
        Constraints = [Constraints, LMI_m <= -1e-6*eye(size(LMI_m))];
    end
    
    % Optimization
    options = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(Constraints, [], options);
    
    if sol.problem == 0
        feasible = true;
        P = value(P_var);
    else
        feasible = false;
        P = NaN(N*n, N*n);
    end
end
