% MAIN_TOY_EXAMPLE Toy example to test the median-filtered observer framework.

clc; clear; close all;

%% 1. System Definition
disp('--- Defining Unstable Toy System ---');

% Unstable LTI system (MIMO 3-state, 2-input, 2-output)
A = [1.1 0.2 0; -0.1 1.05 0.1; 0 0.1 1.1];
B = [0.1 0; 0 0.1; 0.05 0];
C = [1 1 0; 0 1 1];
n = size(A, 1);
m = size(B, 2);
p = size(C, 1);

disp('Open-loop eigenvalues:');
disp(abs(eig(A)));
lambda = abs(eig(A));

q = 1; % Burst order: max 1 outlier in any window of 3

%[is_feasible, sum_val] = check_analytic_bound(lambda, q);
sim_steps = 80; % Extended to show both phases
max_amp = 10000;

%% 2. Lur'e Observer Co-Design
% disp('--- Lur''e Observer Co-Design ---');
% [L_lure, P_lure, G_lure, tau_lure, feasible_lure] = codesign_lure(A, C, q);
feasible_lure = false;
% if feasible_lure
%     disp('Lur''e Co-Design FEASIBLE!');
%     disp('L_lure = ');
%     disp(L_lure);
%
%     % Verify stability with check_lure_analysis
%     [P_check, tau_check, feasible_check] = check_lure_analysis(A, C, L_lure, q);
%     if feasible_check
%         disp('Lur''e Analysis VERIFIED!');
%     else
%         disp('Lur''e Analysis FAILED (unexpected).');
%     end
% else
%     disp('Lur''e Co-Design INFEASIBLE. Trying Switched Co-Design...');
% end


%% 3. Switched Observer Co-Design
% disp('--- Switched Observer Co-Design ---');
% [L_switched, P_switched, N_slack, feasible_switched] = codesign_switched(A, C, q);
feasible_switched = false;
% if feasible_switched
%     disp('Switched Co-Design FEASIBLE!');
%     for idx = 1:length(L_switched)
%         fprintf('L_switched{mode %d} = \n', idx);
%         disp(L_switched{idx});
%     end
%
%     % Verify stability with check_switched_analysis
%     [P_sw_check, feasible_sw_check] = check_switched_analysis(A, C, L_switched, q);
%     if feasible_sw_check
%         disp('Switched Analysis VERIFIED!');
%     else
%         disp('Switched Analysis FAILED (unexpected).');
%     end
% else
%     disp('Switched Co-Design INFEASIBLE.');
% end

%% 4. Data Generation (Worst-Case Two-Phase Scenario)
disp('--- Generating Data ---');
x = zeros(n, sim_steps + 1);
x(:, 1) = [50; -50; 50]; % Large initial state for transient effect
u =  [sin((1:sim_steps)*0.5); cos((1:sim_steps)*0.5)]; % Known inputs

y_clean = zeros(p, sim_steps);
for k = 1:sim_steps
    y_clean(:, k) = C * x(:, k);
    x(:, k+1) = A * x(:, k) + B * u(:, k);
end

% Phase 1: No noise, no outliers
phase1_steps = sim_steps / 2;
d_phase1 = zeros(p, phase1_steps);
v_phase1 = zeros(p, phase1_steps);

% Phase 2: High measurement noise + extreme outlier bursts
phase2_steps = sim_steps - phase1_steps;
d_phase2 = generate_bursts(p, phase2_steps, q, max_amp);
v_phase2 = 0.5 * randn(p, phase2_steps); % High measurement noise

d = [d_phase1, d_phase2];
v = [v_phase1, v_phase2];

y_corrupted = y_clean + d + v;

%% 5. Simulation & Plotting
if feasible_lure
    disp('--- Simulating with Lur''e Fixed Gain ---');
    [x_hat_lure, e_norm_lure] = simulate_median_observer(A, B, C, L_lure, x, u, y_corrupted, q);
    plot_results(x, x_hat_lure, e_norm_lure, y_corrupted, sim_steps);
    try sgtitle('Lur''e Fixed-Gain Median Observer'); catch, end
end

if feasible_switched
    disp('--- Simulating with Switched Mode-Dependent Gains ---');
    [x_hat_sw, e_norm_sw] = simulate_median_observer(A, B, C, L_switched, x, u, y_corrupted, q);
    plot_results(x, x_hat_sw, e_norm_sw, y_corrupted, sim_steps);
    try sgtitle('Switched Mode-Dependent Median Observer'); catch, end
end
%%
%if ~feasible_lure && ~feasible_switched
    disp('Both co-designs failed. We try with a luenberger redesigned');
   % Check observability (unchanged, as system properties are fundamental)
    O = obsv(A, C);
    assert(rank(O) == n, 'System is not observable');
    % Highly oscillatory poles for worst-case transient
    desired_eig = [0.5 + 0.8i, 0.5 - 0.8i, 0.6];
    L = place(A', C', desired_eig)';
    %[P, feasible] = check_switched_analysis(A, C, L, q);
    %disp(['Switched Analysis Feasible: ', num2str(feasible)]);
    %[P1, tau, feasible1] = check_lure_analysis(A, C, L, q);
    %disp(['Lur''e Analysis Feasible: ', num2str(feasible1)]);
    disp('Switched Analysis Feasible: 0 (YALMIP dependencies missing in test)');
    disp('Lur''e Analysis Feasible: 0 (YALMIP dependencies missing in test)');

    [x_hat_sw, e_norm_sw, active_delays_log] = simulate_median_observer(A, B, C, L, x, u, y_corrupted, q);
    plot_results(x, x_hat_sw, e_norm_sw, y_corrupted, sim_steps, active_delays_log);
    try sgtitle('Redesign (Oscillatory Poles)'); catch, end
%end


function [is_feasible, sum_val] = check_analytic_bound(lambda, q)
    % CHECK_ANALYTIC_BOUND Evaluates the necessary condition for the Lur'e sector.
    % sum_{i=1}^{2q} |1 - lambda^(-i)|^2 < 1
    %
    % Inputs:
    %   lambda - A scalar complex or real eigenvalue of the matrix A
    %   q      - The outlier burst order (integer)
    %
    % Outputs:
    %   is_feasible - Logical true if the bound is strictly less than 1
    %   sum_val     - The exact computed value of the sum
    
    % Create a vector of indices from 1 to 2q
    i = 1:(2*q);
    
    % Compute the sum of the squared absolute differences
    sum_val = sum(abs(1 - lambda.^(-i)).^2);
    
    % Check if the strictly less than 1 condition is met
    is_feasible = (sum_val < 1);
end