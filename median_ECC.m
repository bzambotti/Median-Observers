%% Figures of "Outlier Cancellation with Gain-Scheduled Median Filters for Linear Systems"
% Fig. 1 (fig_observers.pdf)
%   (a) estimation error of the Luenberger, stubborn and median observers under
%       bursts of q = 2 consecutive outliers every N = 5 samples;
%   (b) median observer when Definition 1 holds, when it is violated by two
%       bursts of q+1 outliers, and with Gaussian noise added to all samples.
% Fig. 2 (fig_limits.pdf)
%   (a) longest burst q_max certified by the LMIs of Theorem 2 versus T_s
%       (third-order plant), scheduled and single gains;
%   (b) largest certifiable window mu*T_mem versus q for a scalar plant
%       a = exp(mu*Ts), compared with the necessary condition of Proposition 1.
% The numbers quoted in the text are printed in the command window.
% Figures are written to figures/, the LMI data of Fig. 2 to data/fig2_lmi.mat.
% Run from this folder.
%
% Requirements: YALMIP, MOSEK, Control System Toolbox (place),
%               MATLAB R2020a or later (tiledlayout, exportgraphics).

clear; close all; clc;
recompute = false;   % true: recompute the LMI data of Fig. 2 (about 15-30 min with MOSEK)
                     % and save them in data/fig2_lmi.mat; false: load that file
tol = 1e-7;          % LMI margin above which the LMIs are declared feasible

%% Plant (continuous time): eigenvalues {1, -2, -3}
Ac = [0 1 0; 0 0 1; 6 -1 -4];
C  = [1 0 0];
n  = size(Ac, 1);

%% Observer design: Ts = 25 ms, q = 2, lambda = 1
Ts = 0.025; q = 2; lam = 1; T = 6;
N = 2*q + 1;
K = round(T/Ts);
tt = (0:K)*Ts;

[A, Ab, Bb, Cb] = extended_matrices(Ac, C, Ts, q);
[ok, P, Ls] = synthesis(Ab, Bb, Cb, n, exp(-lam*Ts), true);
assert(ok, 'LMIs infeasible');
Llin = place(A', C', exp([-4 -5 -6]*Ts))';      % Luenberger, CT poles {-4,-5,-6}

x0 = [1; 0; 0]; xh0 = zeros(n, 1);
E0 = repmat(x0 - xh0, N, 1);                     % buffer filled with C e_0
kappa = sqrt(max(eig(P))/min(eig(P)));
fprintf('kappa = %.1f, |E0| = %.3f, certified floor 2 kappa |C| |E0| = %.0f\n', ...
        kappa, norm(E0), 2*kappa*norm(C)*norm(E0));
fprintf('|L_i| = %s\n', mat2str(cellfun(@norm, Ls)', 3));

%% Disturbances. Convention: d(t + 2q + 1) = d_t, t = -2q, ..., K
% Nominal: q consecutive outliers every N samples, random sign, amplitude in [200,280]
% (above the certified floor 2 kappa |C| |E0| printed above)
rng(1);
d = zeros(1, K + 1 + 2*q);
for t = 0:K
    if mod(t, N) >= 1 && mod(t, N) <= q
        d(t + 2*q + 1) = (2*(rand > 0.5) - 1) * (200 + 80*rand);
    end
end

% Occasional violation: two bursts of q+1 consecutive outliers with the same sign
dOcc = d;
kViol = zeros(1, 2);
for j = 1:2
    t0 = 80*j;                                   % around t = 2 s and t = 4 s
    k0 = t0 - mod(t0, N) + 1;                    % first sample of a nominal burst
    s0 = sign(dOcc(k0 + 2*q + 1)); if s0 == 0, s0 = 1; end
    dOcc(k0 + 2*q + 1 + (0:q)) = s0*240;         % extend it to q+1 samples
    kViol(j) = k0;
end

% Nominal bursts plus Gaussian measurement noise on all samples
rng(5);
sigmaN = 0.01;
dNoise = d; dNoise(2*q + 1:end) = dNoise(2*q + 1:end) + sigmaN*randn(1, K + 1);

% Persistent violation (text only): q+1 consecutive outliers every N samples
rng(2);
dPer = zeros(1, K + 1 + 2*q);
for t = 0:K
    if mod(t, N) >= 1 && mod(t, N) <= q + 1
        dPer(t + 2*q + 1) = (2*(rand > 0.5) - 1) * (200 + 80*rand);
    end
end

%% Simulations
eL = simulate('lin',  A, C, Ls, Llin, q, d, x0, xh0, K);
eS = simulate('stub', A, C, Ls, Llin, q, d, x0, xh0, K);
[eM,   nbM, bmax] = simulate('med', A, C, Ls, Llin, q, d,    x0, xh0, K);
[eOcc, nbO]       = simulate('med', A, C, Ls, Llin, q, dOcc, x0, xh0, K);
[ePer, nbP]       = simulate('med', A, C, Ls, Llin, q, dPer, x0, xh0, K);
eNoise            = simulate('med', A, C, Ls, Llin, q, dNoise, x0, xh0, K);

fprintf('\nFig. 1(a): final |e|: Luenberger %.2e, stubborn %.2e, median %.2e\n', eL(end), eS(end), eM(end));
fprintf('  median: corrupted selections %d, largest uncorrupted innovation %.2f\n', nbM, bmax);
fprintf('Fig. 1(b): occasional violation: corrupted selections %d, max |e| %.1f, final |e| %.2e\n', ...
        nbO, max(eOcc), eOcc(end));
fprintf('  noise sigma = %.2f: mean |e| over last second %.2e\n', ...
        sigmaN, mean(eNoise(end - round(1/Ts):end)));
fprintf('Persistent violation (text only): corrupted selections %d, mean |e| over last second %.2e\n', ...
        nbP, mean(ePer(end - round(1/Ts):end)));

% Same bursts with amplitude in [50,70], below the certified floor (text only):
% same signs, hence the same selections and the same median error
d50 = d; d50(d ~= 0) = sign(d(d ~= 0)) .* (50 + 20*(abs(d(d ~= 0)) - 200)/80);
[e50, nb50] = simulate('med', A, C, Ls, Llin, q, d50, x0, xh0, K);
fprintf('Amplitude [50,70]: corrupted selections %d, final |e| %.2e\n', nb50, e50(end));

% Stubborn observer: residual error versus Theta (text only)
for Theta = [0.01 0.1 0.5 2]
    eTh = simulate('stub', A, C, Ls, Llin, q, d, x0, xh0, K, Theta);
    fprintf('Stubborn, Theta = %.2f: mean |e| over last second %.2f\n', Theta, mean(eTh(end - round(1/Ts):end)));
end

% Isolated outliers separated by two clean samples, median designed with q = 1 (text only)
q1 = 1; lam1 = 3;
[~, Ab1, Bb1, Cb1] = extended_matrices(Ac, C, Ts, q1);
[ok1, P1, Ls1] = synthesis(Ab1, Bb1, Cb1, n, exp(-lam1*Ts), true);
assert(ok1, 'LMIs infeasible for q = 1');
kappa1 = sqrt(max(eig(P1))/min(eig(P1)));
rng(1);
d1 = zeros(1, K + 1 + 2*q);                     % same indexing as for q = 2
for t = 0:K
    if mod(t, 3) == 1
        d1(t + 2*q + 1) = (2*(rand > 0.5) - 1) * (200 + 80*rand);
    end
end
d1q = d1(2*q - 2*q1 + 1:end);                   % re-index for q1: d1q(t + 2*q1 + 1) = d_t
e1 = simulate('med', A, C, Ls1, Llin, q1, d1q, x0, xh0, K);
fprintf('Isolated outliers: floor (q = 1) %.0f, median final |e| %.1e\n', 2*kappa1*norm(C)*sqrt(3), e1(end));
for Theta = [0.02 0.1 0.2]
    eTh = simulate('stub', A, C, Ls, Llin, q, d1, x0, xh0, K, Theta);
    fprintf('  stubborn, Theta = %.2f: mean |e| over last second %.2e\n', Theta, mean(eTh(end - round(1/Ts):end)));
end

% No outliers: price of robustness (text only)
d0 = zeros(1, K + 1 + 2*q);
fprintf('No outliers: final |e| median %.1e, Luenberger %.1e\n', ...
        simulate_last('med', A, C, Ls, Llin, q, d0, x0, xh0, K), ...
        simulate_last('lin', A, C, Ls, Llin, q, d0, x0, xh0, K));

% Small outliers, below the magnitude condition (text only)
rng(1);
dSm = zeros(1, K + 1 + 2*q);
for t = 0:K
    if mod(t, N) >= 1 && mod(t, N) <= q
        dSm(t + 2*q + 1) = (2*(rand > 0.5) - 1) * (0.05 + 0.15*rand);
    end
end
[eSm, nbSm] = simulate('med', A, C, Ls, Llin, q, dSm, x0, xh0, K);
fprintf('Small outliers in [0.05,0.2]: corrupted selections %d, final |e| %.2e (nominal %.2e)\n', ...
        nbSm, eSm(end), eM(end));

% Measurement noise on all samples (text only)
for sigma = 0.1                                  % sigma = 0.01 is eNoise above
    rng(5);
    dN = d; dN(2*q + 1:end) = dN(2*q + 1:end) + sigma*randn(1, K + 1);
    eN = simulate('med', A, C, Ls, Llin, q, dN, x0, xh0, K);
    fprintf('Noise sigma = %.2f: mean |e| of the median observer over last second %.2e\n', ...
            sigma, mean(eN(end - round(1/Ts):end)));
end

%% Fig. 2(a): q_max versus Ts (third-order plant)
TsList = [0.0125 0.025 0.05 0.1 0.15];
qList = 1:6;
dataFile = fullfile('data', 'fig2_lmi.mat');
if recompute
    qS = zeros(size(TsList)); qU = qS;           % q_max, scheduled / single gain
    for j = 1:numel(TsList)
        qS(j) = q_max(Ac, C, TsList(j), true,  tol);
        qU(j) = q_max(Ac, C, TsList(j), false, tol);
    end
    aS = zeros(size(qList)); aU = aS;            % a_max, scheduled / single gain
    for j = 1:numel(qList)
        aS(j) = scalar_ceiling(qList(j), true,  tol);
        aU(j) = scalar_ceiling(qList(j), false, tol);
        fprintf('  q = %d: a_max scheduled = %.4f, single = %.4f\n', qList(j), aS(j), aU(j));
    end
    if ~exist('data', 'dir'), mkdir('data'); end
    save(dataFile, 'TsList', 'qList', 'qS', 'qU', 'aS', 'aU', 'tol');
elseif exist(dataFile, 'file')
    S = load(dataFile);
    assert(isequal(S.TsList, TsList) && isequal(S.qList, qList), ...
           'data/fig2_lmi.mat was computed for other Ts or q: set recompute = true');
    qS = S.qS; qU = S.qU; aS = S.aS; aU = S.aU;
else
    % Values from the last run before data/fig2_lmi.mat existed (Oct 2026).
    % Run once with recompute = true to replace them with saved results.
    warning('data/fig2_lmi.mat not found: using stored values');
    qS = [4 3 2 1 0];               qU = [3 2 1 1 0];
    aS = [1.1547 1.0515 1.0257 1.0154 1.0103 1.0073];
    aU = [1.1382 1.0430 1.0207 1.0121 1.0080 1.0056];
end
fprintf('\nFig. 2(a):\n');
for j = 1:numel(TsList)
    fprintf('  Ts = %5.1f ms: q_max scheduled = %d (burst %.0f ms), single = %d\n', ...
            1000*TsList(j), qS(j), 1000*qS(j)*TsList(j), qU(j));
end

%% Fig. 2(b): scalar plant a = exp(mu*Ts), largest mu*T_mem = 2q ln(a_max)
nec = 2*qList.*log(1 + 1./(2*qList));           % Proposition 1
muS = 2*qList.*log(aS);
muU = 2*qList.*log(aU);
fprintf('\nFig. 2(b): mu*T_mem\n');
fprintf('  necessary: %s\n  scheduled: %s\n  single:    %s\n', ...
        mat2str(nec, 3), mat2str(muS, 3), mat2str(muU, 3));

%% Style
BLUE = hex2rgb('#2a78d6'); LBLUE = hex2rgb('#86b6ef'); DBLUE = hex2rgb('#0d366b');
ORANGE = hex2rgb('#eb6834'); AQUA = hex2rgb('#1baf7a');
GRAY = hex2rgb('#52514e'); BAND = hex2rgb('#c9c8c3');

%% Fig. 1
fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 20, 9]);
tiledlayout(fig1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile; hold(ax1, 'on'); set(ax1, 'YScale', 'log');
h1 = plot(ax1, tt, eL, '--', 'Color', ORANGE, 'LineWidth', 2);
h2 = plot(ax1, tt, eS, '-.', 'Color', AQUA,   'LineWidth', 2);
h3 = plot(ax1, tt, eM, '-',  'Color', BLUE,   'LineWidth', 2);
legend(ax1, [h1 h2 h3], {'Luenberger', 'Stubborn', 'Median'}, ...
       'Location', 'southwest', 'FontSize', 14, 'Box', 'off');
ylabel(ax1, '$|e|$', 'Interpreter', 'latex');


ax2 = nexttile; hold(ax2, 'on'); set(ax2, 'YScale', 'log');
g3 = plot(ax2, tt, eNoise, '-.', 'Color', LBLUE, 'LineWidth', 2);
g2 = plot(ax2, tt, eOcc,   '--', 'Color', DBLUE, 'LineWidth', 2);
g1 = plot(ax2, tt, eM,     '-',  'Color', BLUE,  'LineWidth', 2);
lgd = legend(ax2, [g1 g2 g3], {'Def. 1 holds', 'Def. 1 violated', 'Outliers + Persistent noise'}, ...
       'Location', 'southwest', 'FontSize', 14, 'Box', 'off');
lgd.ItemTokenSize = [40, 18];                    % longer line samples


for ax = [ax1 ax2]
    xlim(ax, [0 T]); ylim(ax, [1e-6 1e3]); yticks(ax, [1e-6 1e-3 1 1e3]);
    xlabel(ax, '$t$ [s]', 'Interpreter', 'latex');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontName', 'Times', 'FontSize', 16, 'TickLabelInterpreter', 'latex');
end
yticklabels(ax2, {});                            % shared y-axis with (a)
if ~exist('figures', 'dir'), mkdir('figures'); end
exportgraphics(fig1, fullfile('figures', 'fig_observers.pdf'), 'ContentType', 'vector');

%% Fig. 2
fig2 = figure('Units', 'centimeters', 'Position', [2, 12, 20, 9]);
tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

bx1 = nexttile; hold(bx1, 'on');
k2 = plot(bx1, 1000*TsList, qU, '--s', 'Color', BLUE, 'LineWidth', 2, ...
          'MarkerSize', 10, 'MarkerFaceColor', 'w');
k1 = plot(bx1, 1000*TsList, qS, '-o', 'Color', BLUE, 'LineWidth', 2, ...
          'MarkerSize', 8, 'MarkerFaceColor', BLUE);
legend(bx1, [k1 k2], {'scheduled', 'single gain'}, ...
       'Location', 'northeast', 'FontSize', 14, 'Box', 'off');
xlim(bx1, [0 160]); ylim(bx1, [-0.3 max(qS) + 0.6]); yticks(bx1, 0:max(qS));
xlabel(bx1, '$T_s$ [ms]', 'Interpreter', 'latex');
ylabel(bx1, '$q_{\max}$', 'Interpreter', 'latex');

bx2 = nexttile; hold(bx2, 'on');
m0 = plot(bx2, qList, nec, ':^', 'Color', GRAY, 'LineWidth', 2, ...
          'MarkerSize', 8, 'MarkerFaceColor', GRAY);
m2 = plot(bx2, qList, muU, '--s', 'Color', BLUE, 'LineWidth', 2, ...
          'MarkerSize', 10, 'MarkerFaceColor', 'w');
m1 = plot(bx2, qList, muS, '-o', 'Color', BLUE, 'LineWidth', 2, ...
          'MarkerSize', 8, 'MarkerFaceColor', BLUE);
legend(bx2, [m0 m1 m2], {'necessary', 'LMI, scheduled', 'LMI, single gain'}, ...
       'Location', 'east', 'FontSize', 14, 'Box', 'off');
xlim(bx2, [0.7 max(qList) + 0.3]); ylim(bx2, [0 1.05]); xticks(bx2, qList);
xlabel(bx2, '$q$', 'Interpreter', 'latex');
ylabel(bx2, '$\mu T_{\mathrm{mem}}$', 'Interpreter', 'latex');


for ax = [bx1 bx2]
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontName', 'Times', 'FontSize', 16, 'TickLabelInterpreter', 'latex');
end
exportgraphics(fig2, fullfile('figures', 'fig_limits.pdf'), 'ContentType', 'vector');

%% ========================================================================
%  Local functions
%  ========================================================================
function [Ab, Bb, Cb] = extended_from_A(A, C, q)
% Extended matrices Abar, Bbar, Cbar_i (cell index i+1, i = 0..2q) for a
% discrete-time pair (A, C).
n = size(A, 1); N = 2*q + 1;
Ab = zeros(N*n);
Ab(1:n, 1:n) = A;
for j = 2:N                                      % shift of the past errors
    Ab((j-1)*n + (1:n), (j-2)*n + (1:n)) = eye(n);
end
Bb = [eye(n); zeros((N-1)*n, n)];
Cb = cell(N, 1);
for i = 0:2*q                                    % Cbar_i E_k = C e_{k-i}
    Ci = zeros(1, N*n);
    Ci(i*n + (1:n)) = C;
    Cb{i+1} = Ci;
end
end

function [A, Ab, Bb, Cb] = extended_matrices(Ac, C, Ts, q)
% Exact discretization A = expm(Ac*Ts) and extended matrices
A = expm(Ac*Ts);
[Ab, Bb, Cb] = extended_from_A(A, C, q);
end

function [ok, Pv, Ls] = synthesis(Ab, Bb, Cb, n, alpha, scheduled)
% LMIs of Theorem 2 with slack G = [G1 G2; 0 G3] and Y_i = G1 L_i.
% scheduled = false imposes Y_i = Y (single gain).
yalmip('clear');
M = size(Ab, 1); N = numel(Cb); ep = 1e-7;
P  = sdpvar(M, M);
G1 = sdpvar(n, n, 'full'); G2 = sdpvar(n, M - n, 'full'); G3 = sdpvar(M - n, M - n, 'full');
G  = [G1, G2; zeros(M - n, n), G3];
Y = cell(N, 1);
if scheduled
    for i = 1:N, Y{i} = sdpvar(n, 1); end
else
    Ys = sdpvar(n, 1);
    for i = 1:N, Y{i} = Ys; end
end
cons = [P >= ep*eye(M)];
for i = 1:N
    X = G*Ab - Bb*Y{i}*Cb{i};
    cons = [cons, [-alpha^2*P, X'; X, P - G - G'] <= -ep*eye(2*M)]; %#ok<AGROW>
end
sol = optimize(cons, [], sdpsettings('solver', 'mosek', 'verbose', 0));
ok = false; Pv = []; Ls = cell(N, 1);
if sol.problem ~= 0, return; end
Pv = value(P); Pv = (Pv + Pv')/2;
G1v = value(G1);
for i = 1:N
    Ls{i} = G1v \ value(Y{i});
end
ok = true;
end

function t = lmi_margin(Ab, Bb, Cb, n, scheduled)
% Largest t such that the LMIs of Theorem 2 with alpha = 1 hold with margin t,
% under the normalization trace(P) = dim(P). The LMIs are feasible for some
% alpha < 1 if and only if t > 0. Returns NaN if the solver fails.
yalmip('clear');
M = size(Ab, 1); N = numel(Cb);
P  = sdpvar(M, M);
G1 = sdpvar(n, n, 'full'); G2 = sdpvar(n, M - n, 'full'); G3 = sdpvar(M - n, M - n, 'full');
G  = [G1, G2; zeros(M - n, n), G3];
Y = cell(N, 1);
if scheduled
    for i = 1:N, Y{i} = sdpvar(n, 1); end
else
    Ys = sdpvar(n, 1);
    for i = 1:N, Y{i} = Ys; end
end
tv = sdpvar(1);
cons = [P >= 0, trace(P) == M, tv <= 1];
for i = 1:N
    X = G*Ab - Bb*Y{i}*Cb{i};
    cons = [cons, [-P, X'; X, P - G - G'] <= -tv*eye(2*M)]; %#ok<AGROW>
end
sol = optimize(cons, -tv, sdpsettings('solver', 'mosek', 'verbose', 0));
if sol.problem == 0
    t = value(tv);
else
    t = NaN;
end
end

function qm = q_max(Ac, C, Ts, scheduled, tol)
% Largest q for which the LMIs of Theorem 2 are feasible for some alpha < 1
n = size(Ac, 1);
qm = 0;
while qm < 12
    [~, Ab, Bb, Cb] = extended_matrices(Ac, C, Ts, qm + 1);
    t = lmi_margin(Ab, Bb, Cb, n, scheduled);
    if ~(t > tol), break; end
    qm = qm + 1;
end
end

function amax = scalar_ceiling(q, scheduled, tol)
% Largest a > 1 for which the LMIs of Theorem 2 are feasible for the scalar
% plant A = a, C = 1 (bisection; a < 1 + 1/(2q) by Proposition 1)
lo = 1; hi = 1 + 1/(2*q);
for it = 1:22
    mid = (lo + hi)/2;
    [Ab, Bb, Cb] = extended_from_A(mid, 1, q);
    if lmi_margin(Ab, Bb, Cb, 1, scheduled) > tol
        lo = mid;
    else
        hi = mid;
    end
end
amax = lo;
end

function [err, nbad, bmax] = simulate(obs, A, C, Ls, Llin, q, d, x0, xh0, K, Theta, tau)
% obs: 'lin' (Luenberger), 'stub' (finite-memory stubborn), 'med' (median).
% d(t + 2q + 1) = d_t. Returns |e(t_k)|, k = 0..K, the number of corrupted
% samples selected by the median, and the largest uncorrupted innovation |C e_k|.
if nargin < 11, Theta = 0.1; end
if nargin < 12, tau = 3; end
N = 2*q + 1;
x = x0; xh = xh0;
innov0 = C*(x0 - xh0);
buf  = innov0*ones(1, N - 1);                    % past innovations (median)
sbuf = abs(innov0)*ones(1, tau);                 % past |innovations| (stubborn)
err = zeros(1, K + 1); nbad = 0; bmax = abs(innov0);
for k = 0:K
    ce = C*(x - xh);                             % uncorrupted innovation C e_k
    yt = ce + d(k + 2*q + 1);                    % innovation y~_k
    err(k + 1) = norm(x - xh);
    bmax = max(bmax, abs(ce));
    switch obs
        case 'lin'
            inj = Llin*yt;
        case 'stub'
            Delta = Theta*sum(sbuf);             % threshold from past innovations
            inj = Llin*min(max(yt, -Delta), Delta);
            sbuf = [abs(yt), sbuf(1:end-1)];
        case 'med'
            s  = [yt, buf];                      % sigma_k = (y~_k, ..., y~_{k-2q})
            ss = sort(s);
            med = ss(q + 1);
            m = find(s == med, 1) - 1;           % smallest delay attaining the median
            if d(k - m + 2*q + 1) ~= 0, nbad = nbad + 1; end
            inj = Ls{m + 1}*med;
            buf = [yt, buf(1:end-1)];
    end
    xh = A*xh + inj;
    x  = A*x;
end
end

function v = simulate_last(varargin)
% Final value of |e| returned by simulate
e = simulate(varargin{:});
v = e(end);
end

function c = hex2rgb(h)
c = sscanf(h(2:end), '%2x%2x%2x', [1 3])/255;
end