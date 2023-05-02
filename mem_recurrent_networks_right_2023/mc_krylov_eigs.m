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
method_A = 'normal';
rho = 0.9;

% Method to sample the input mask
method_C = 'unit';

switch method_A
    case 'normal'
        A_ = randn(N);
    case 'uniform'
        A_ = (rand(N)-.5) * 2;
    case 'sparse_normal'
        A_ = full(sprandn(N,N,0.1));
    case 'sparse_normal_cond'
        A_ = full(sprandn(N,N,0.1,0.8));
    case 'ortho'
        A_ = randortho(N);
    case 'cyclic'
        A_ = eye(N);
        A_ = A_(:,[N, 1:N-1]);
    otherwise
        error("No appropriate sampling rule for A chosen.")
end

A = A_ / max(abs(eig(A_))) * rho;
% A = A_ / max(svd(A_)) * rho;

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

A_abseigs = [1; sort(abs(eig(A)), 'desc')];
Powiter_fac = A_abseigs(2) / A_abseigs(1);
A_weakvals =  A_abseigs(A_abseigs < rho);

%% Simulation

X = zeros(N,N*5);
Ac = eye(N)*C;
for i = 1:size(X,2)
    X(:,i) = Ac;
    Ac = A * Ac;
end
[MCs_arn, Q, ~] = mc_arnoldi(A, C);

d_Acv     = zeros(1,round(1.5*N));
d_Acv_arn = zeros(1,round(1.5*N));
d_Acv_ort = zeros(1,round(1.5*N));
d_eigs    = zeros(1,round(1.5*N));
for i = 2:size(d_Acv,2) 
    % compute difference from projection
    Hi = X(:,1:i-1);
    pHi = Hi * pinv(Hi);
    % pHi = Hi / (Hi' * Hi + 10^-10*eye(i-1)) * Hi';
    
    Qi = Q(:,1:min(i-1,N));
    pHi_arn = Qi * pinv(Qi);
    
    [Ui, Di, Wi] = svd(Hi,'econ');
    pHi_svd = Ui * Ui';
    
    vi = X(:,i); % / norm(X(:,i));
    vHi = pHi * vi;
    d_Acv(i)     = norm(vi - pHi * vi);
    d_Acv_arn(i) = norm(vi - pHi_arn * vi);
    d_Acv_ort(i) = norm(vi - pHi_svd * vi);
    
    if i < N+1
        d_eigs(i) = prod(A_abseigs(1:i));
    else
        d_eigs(i) = prod(A_abseigs);
    end
end

%% Plot

fig = figure(101);
tiledlayout(1, 1, 'TileSpacing', 'none')
nexttile
% plot(0:size(d_Acv,2)-1, d_Ac_v, 'Marker', '.');
% p1 = semilogy(0:size(d_Acv,2)-1, d_Acv, 'Marker', '.', 'LineWidth', 1);
% hold on
p2 = semilogy(0:size(d_Acv,2)-1, d_Acv_arn, ...
                    'Marker', '.', 'LineWidth', 1);
hold on
p3 = semilogy(0:size(d_Acv,2)-1, d_Acv_ort, ...
                    'Marker', '.', 'LineWidth', 1);
p4 = semilogy(0:size(d_Acv,2)-1, d_eigs, 'k--');
semilogy(0:size(d_Acv,2)-1, (rho).^(0:size(d_Acv,2)-1), 'k:')
hold off
xline(N, 'Label', "N = " + N, ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right')
xline(rank(X), 'Label', "rank(K_N) = " + rank(X), 'Color', '#7E2F8E', ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left')
yline(eps, 'k-', "eps", 'LabelHorizontalAlignment', 'left')
xlim([0, size(d_Acv,2)-1])
ylim([10^-25 1])
xlabel('$j$', 'interpreter', 'latex')
ylabel('$||\boldmath{\theta}_j||$', 'interpreter', 'latex')
grid
% title("Krylov subspace squeeze")
legend([p2, p3, p4], "Arnoldi", "Ortho", "Eigen", ...
        'Location', 'northeast', 'Orientation', 'vertical')

% Save figure
if SAVE_PLOT_FLAG
printpdf(fig, ...
    join(["figures/", "plot_subs_krylov_eigs_", ...
            "_A=", string(method_A), ...
            "_rho=", string(rho), ...
            "_C=", string(method_C), ".pdf"],""), ...
    [0, 0, 12, 10])
    disp("Plot saved!")
end
   
% #####