function [L, P, G, tau, feasible] = codesign_lure(A, C, q)
    % CODESIGN_LURE Synthesizes a robust median-filtered observer using the Lur'e framework.
    % 
    % Inputs:
    %   A - System matrix (n x n)
    %   C - Output matrix (p x n)
    %   q - Outlier burst order (scalar)
    %
    % Outputs:
    %   L - Synthesized observer gain (n x p)
    %   P - Lyapunov matrix (Nn x Nn)
    %   G - Slack matrix
    %   tau - S-procedure scalar
    %   feasible - Boolean indicating if the LMI was solved successfully
    
    n = size(A, 1);
    p = size(C, 1);
    N = 2 * q + 1;
    
    % Define \bar{A}, \bar{B}, \bar{C}
    A_bar = zeros(N*n, N*n);
    A_bar(1:n, 1:n) = A;
    for i = 1:(N-1)
        A_bar(i*n+1:(i+1)*n, (i-1)*n+1:i*n) = eye(n);
    end
    
    B_bar = zeros(N*n, n);
    B_bar(1:n, 1:n) = eye(n);
    
    C_bar = zeros(p, N*n);
    C_bar(1:p, 1:n) = C;
    
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
    G1 = sdpvar(n, n, 'full');
    G2 = sdpvar(n, 2*q*n, 'full');
    G3 = sdpvar(2*q*n, 2*q*n, 'full');
    G = [G1, G2; zeros(2*q*n, n), G3];
    Y = sdpvar(n, p, 'full');
    tau = sdpvar(1, 1);
    
    % LMI definition
    M11 = -P + tau * Q_o;
    M13 = A_bar' * G' - C_bar' * Y' * B_bar';
    M22 = -tau * eye(p);
    M23 = Y' * B_bar';
    M33 = P - G - G';
    
    LMI = [M11,           zeros(N*n, p), M13;
           zeros(p, N*n), M22,           M23;
           M13',          M23',          M33];
           
    Constraints = [LMI <= -1e-6*eye(size(LMI)), P >= 1e-6*eye(N*n), tau >= 1e-6];
    
    % Optimization
    options = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(Constraints, [], options);
    
    if sol.problem == 0
        feasible = true;
        L = value(G1) \ value(Y);
        P = value(P);
        G = value(G);
        tau = value(tau);
    else
        feasible = false;
        L = NaN(n, p);
        P = NaN(N*n, N*n);
        G = NaN(N*n, N*n);
        tau = NaN;
    end
end
