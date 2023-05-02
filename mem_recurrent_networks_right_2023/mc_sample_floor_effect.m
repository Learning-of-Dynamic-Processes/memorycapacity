%
%
%

clear
rng('default')
 
% Save plots?
SAVE_PLOT_FLAG = true;

%% Parameters

N = 100;            % ESN size
M = floor(5*N);     % Number of MCs to compute

% Method to sample the connectivity matrix, spectral radius
method_A = 'ortho';
rho = 0.9;

% Method to sample the input mask
method_C = 'normal';

% Sample size to simulate
T = 1000:500:10000;


%% ESN model

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

sims = struct([]);

discard = 100;

for j = 1:length(T)
    S = T(j);
    
    % State linear
    E = randn(1, S+discard);
    X_lin = zeros(N, S+discard);
    X_lin(:,1) = C*E(:,1);
    for t = 1:S+discard-1
        X_lin(:,t+1) = A*X_lin(:,t) + C*E(:,t);
    end
    X_lin = X_lin(:, 1+discard:end);
    z = E(discard:end-1);
    
    Gx = X_lin*X_lin';
    g  = z*z';
    
    MC_tau = nan(1, M);
    % Compute MC
    for i = 1:M
        MCov_i = sum(X_lin(:,i:end) .* z(1:end-i+1), 2);
        % MC_tau(i) = (MCov_i' / (Gx + 1e-9*eye(N)) * MCov_i) / g;
        % MC_tau(i) = (MCov_i' * pinv(Gx) * MCov_i) / g;
        MC_tau(i) = (MCov_i' / (Gx) * MCov_i) / g;
    end
    
    sims{j}.MC_tau = MC_tau;
    sims{j}.MC     = sum(MC_tau) / N;
end

% Unwrap MC
MCs = zeros(1, length(T));
for j = 1:length(T)
    MCs(j) = sims{j}.MC;
end

% Subspace MC
X = zeros(N,M);
Ac = eye(N)*C;
for i = 1:size(X,2)
    X(:,i) = Ac;
    Ac = A * Ac;
end
[U, D, V] = svd(X,'econ');
sMC = diag(V * V');

% Difference of MC_tau
Df_MC_tau = zeros(1, length(T));
for j = 1:length(T)
    Df_MC_tau(j) = norm(sims{j}.MC_tau(:) - sMC(:));
end
   
%% Plot

% cmap = flip(winter(length(T)), 1);
cmap = flip(turbo(length(T)+4), 1);
cmap = cmap(3:end-2,:);

fig = figure();

t = tiledlayout(1, 2, 'Padding', 'compact');
% title(t, sprintf("Inconsistent sample MC estimation - N = %d", N))

nexttile
colormap(cmap)
% figure
hold on
for j = 1:length(T)
    plot(1:round(2*N), sims{j}.MC_tau(1:round(2*N)), 'Color', cmap(j,:))
end
% plot(1:round(2*N), sMC(1:round(2*N)), 'Color', '#D95319')
grid on
ylim([0, 1.05])
yline(1, 'k--')
xline(N, 'Label', "N = " + N, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right')
xlabel("$\tau$", 'interpreter', 'latex')
ylabel("$\widehat{\textnormal{MC}}_\tau(T)$", 'interpreter', 'latex')
% caxis([min(T), max(T)])
% colorbar('Ticks', [min(T), mean(T), max(T)], ...
%          'TickLabels', {num2str(min(T)), num2str(mean(T)), num2str(max(T))}, ...
%          'Location', 'southoutside');
     
nexttile
hold on
b = bar(T, MCs);
b.FaceColor = 'flat';
for j = 1:length(T)
    b.CData(j,:) = cmap(j,:);
end
grid on
yline(1, 'k--')
xlabel("$T$", 'interpreter', 'latex')
ylabel("$\widehat{\textnormal{MC}}(T)$", 'interpreter', 'latex')
% cb.Layout.Tile = 'east';

% Save figure
if SAVE_PLOT_FLAG
printpdf(fig, ...
    join(["figures/", "plot_sample_floor", ...
            "_N=", string(N), ...
            "_A=", string(method_A), ...
            "_C=", string(method_C), ".pdf"],""), ...
    [0, 0, 24, 9])
end
    
% #####