function [P, tau, feasible] = check_lure_analysis(A, C, L, q)
    % CHECK_LURE_ANALYSIS Verifies the stability of a given observer gain using the Lur'e framework.
    % 
    % Inputs:
    %   A - System matrix (n x n)
    %   C - Output matrix (p x n)
    %   L - Observer gain (n x p)
    %   q - Outlier burst order (scalar)
    %
    % Outputs:
    %   P - Lyapunov matrix (Nn x Nn)
    %   tau - S-procedure scalar
    %   feasible - Boolean indicating if the LMI was solved successfully
    
    n = size(A, 1);
    p = size(C, 1);
    N = 2 * q + 1;
    
    % Define \mathcal{A}_o, \mathcal{B}_o
    A_o = zeros(N*n, N*n);
    A_o(1:n, 1:n) = A - L*C;
    for i = 1:(N-1)
        A_o(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = eye(n);
    end
    
    B_o = zeros(N*n, p);
    B_o(1:n, 1:p) = L;
    
    % Define \mathcal{Q}_o
    Q_o = zeros(N*n, N*n);
    CTC = C' * C;
    Q_o(1:n, 1:n) = 2 * q * CTC;
    for i = 1:(2*q)
        idx = i*n + 1 : (i+1)*n;
        Q_o(1:n, idx) = -CTC;
        Q_o(idx, 1:n) = -CTC;
        Q_o(idx, idx) = CTC;
    end
    
    % YALMIP variables
    P = sdpvar(N*n, N*n, 'symmetric');
    tau = sdpvar(1, 1);
    
    % LMI definition
    M11 = A_o' * P * A_o - P + tau * Q_o;
    M12 = A_o' * P * B_o;
    M22 = B_o' * P * B_o - tau * eye(p);
    
    LMI = [M11, M12;
           M12', M22];
           
    Constraints = [LMI <= -1e-6*eye(size(LMI)), P >= 1e-6*eye(N*n), tau >= 1e-6];
    
    % Optimization
    options = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(Constraints, [], options);
    
    if sol.problem == 0
        feasible = true;
        P = value(P);
        tau = value(tau);
    else
        feasible = false;
        P = NaN(N*n, N*n);
        tau = NaN;
    end
end
