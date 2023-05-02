%
%
%

clear
rng('default')
 
% Save plots?
SAVE_PLOT_FLAG = true;

%% Parameters

N = 100;

% Method to sample the connectivity matrix, spectral radius
method_A = 'sparse_normal';
rho = 0.9;

% Method to sample the input mask
method_C = 'sparse_normal';

switch method_A
    case 'normal'
        A_ = randn(N);
    case 'uniform'
        A_ = (rand(N)-.5) * 2;
    case 'sparse_normal'
        A_ = full(sprandn(N,N,0.1));
    case 'sparse_normal_cond'
        A_ = full(sprandn(N,N,0.1,0.7));
    case 'ortho'
        A_ = randortho(N);
    case 'cyclic'
        A_ = eye(N);
        A_ = A_(:,[N, 1:N-1]);
    case 'delay'
        A_ = toeplitz([0, 1, zeros(1,N-2)], zeros(1,N));
    otherwise
        error("No appropriate sampling rule for A chosen.")
end

maxEigA_ = max(abs(eig(A_)));
if maxEigA_ > 0
    A = A_ / max(abs(eig(A_))) * rho;
    % A = A_ / max(svd(A_)) * rho;
else
    A = A_ * rho;
end

switch method_C
    case 'normal'
        C = randn(N,1);
    case 'sparse_normal'
        C = full(sprandn(N,1,0.1));
    case 'uniform'
        C = (rand(N,1)-.5)*2;
    case 'unit'
        C = ones(N,1);
    case 'base'
        C = [1; zeros(N-1,1)];
    otherwise
        error("No appropriate sampling rule for C chosen.")
end

C = C / norm(C); % normalize input mask

%% Simulation

X = zeros(N,N*5);
Ac = eye(N)*C;
for i = 1:size(X,2)
    X(:,i) = Ac;
    Ac = A * Ac;
end

% Naive MC computation
MCj_naive = diag(X' * pinv(X'));
MC = sum(MCj_naive) / N;

% Replications for OSM+
B = 1000;
MCj_subsp = zeros(N*2, B);
for b = 1:B
    % Change rnd seed
    rng(1039493+b*17)
    
    switch method_C
        case 'normal'
            C = randn(N,1);
        case 'uniform'
            C = (rand(N,1)-.5)*2;
        case 'sparse_normal'
            C = full(sprandn(N,1,0.1));
        case 'sparse_uniform'
            C = full(sprand(N,1,0.1));
        case 'unit'
            C = ones(N,1);
        case 'base'
            C = [1; zeros(N-1,1)];
        otherwise
            error("No appropriate sampling rule for C chosen.")
    end
    
    C = C / norm(C);    % normalize input mask
    
    X = zeros(N,N*5);
    Ac = eye(N)*C;
    for i = 1:size(X,2)
        X(:,i) = Ac;
        Ac = A * Ac;
    end
    % Subspace MC
    [U, D, V] = svd(X,'econ');
    W = V;
    % W = qr(X);
    % W = gramschmidt(X');
    sMCi = diag(W * W');
    sMC = sum(sMCi) / N;
    MCj_subsp(:,b) = sMCi(1:N*2);
    
end
MCj_subsp_m = mean(MCj_subsp, 2);
    
%% Plot

fig = figure();
tiledlayout(1, 1, 'Padding', 'Compact');

nexttile
p1 = plot(1:round(1.5*N), MCj_naive(1:round(1.5*N)), 'Marker', '.');
hold on
p2 = plot(1:round(1.5*N), MCj_subsp_m(1:round(1.5*N)), 'Marker', '.');
hold off
xline(N, 'Label', "N = " + N, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right')
grid on
grid minor
ylim([0 1.05])
xlabel("$\tau$", 'interpreter', 'latex')
ylabel("$\widehat{\textnormal{MC}}_\tau$", 'interpreter', 'latex')
title("Memory Capacity Curves")
legend([p1, p2], "naive", "OSM+", ...
        'Location', 'southoutside', 'Orientation', 'horizontal')
subtitle( sprintf("Total MC:  naive = %g,  OSM/OSM+ = %g", ...
    100 * MC, 100 * sMC))

% Save figure
if SAVE_PLOT_FLAG
printpdf(fig, ...
    join(["figures/", "plot_mc_", ...
            "_N=", string(N), ...
            "_A=", string(method_A), ...
            "_C=", string(method_C), ".pdf"],""), ...
    [0, 0, 12, 10])
end
    
% #####