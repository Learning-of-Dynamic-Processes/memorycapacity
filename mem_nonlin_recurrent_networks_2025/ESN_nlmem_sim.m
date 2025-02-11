%
%
%

clear
rng('default')

T = 1e4;
N = 30;

fun = @tanh;

% FLAGS
SAVE_FIGURE = false;


%% Sampling options

% Method to sample the connectivity matrix, spectral radius
method_A = 'ortho';
rho = 0.95;

% Method to sample the input mask
method_C = 'normal';

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
    otherwise
        error("No appropriate sampling rule for A chosen.")
end

%A = A_ / max(abs(eig(A_))) * rho;
A = A_ / max(svd(A_)) * rho;

switch method_C
    case 'normal'
        C_ = randn(N,1);
    case 'uniform'
        C_ = (rand(N,1)-.5)*2;
    case 'unit'
        C_ = ones(N,1);
    case 'base'
            C_ = [1; zeros(N-1,1)];
    otherwise
        error("No appropriate sampling rule for C chosen.")
end

%% Simulate
Rad = randi(2, T, 1) * 2 - 3;
% Rad = randn(T, 1);

lambda = 0; % regularization
% true bound
gamma_grid_min = 0.4*(1-rho)/(max(abs(C_)));
% bound for simulations
% gamma_grid_min = 0.4*(1-rho)/(sqrt(N)*max(abs(C_)));
D = 4;
% looser bound used for simulations
gamma_grid_max = (D+max(abs(A(:)))*N)/min(abs(C_));

gamma_grid = logspace(log10(gamma_grid_min),log10(gamma_grid_max),17);
% gamma_grid = [10^(-3), 10^floor(log10(gamma_grid_min)), gamma_grid];

L = length(gamma_grid);
MC_by_gamma = zeros(L,1);
for i = 1:L
    % NOTE: instead of rescaling the Rademacher r.v.s, for simplicity
    %       we rescale the input matrix. This is equivalent.
    % C = C_ / abs(max(C_)) * gamma_grid(i);
    C = C_ * gamma_grid(i);
    % Collect states
    discard = 100;
    X = zeros(N, T);
    X(:,1) = C*Rad(1);
    for t = 2:T
        X(:,t) = fun(A*X(:,t-1) + C*Rad(t));
    end
    % Rescale nonlinear states (without loss of generality)
    X = X(:,(discard+1):end);
    X = X/min(gamma_grid(i), 1);  
    Y = Rad((discard+1):end)';

    MC = zeros(2*N,1);
    for j = 1:(2*N)
        k = j-1;
        X_j = X(:,k+1:end);
        Y_j = Y(:,1:end-k);
        
        % Regression
        if lambda > 0
            W0_j = (X_j*X_j' + T*lambda*eye(K)) \ (X_j*Y_j');
            a_j = mean(Y_j,2) - W0_j' * mean(X_j,2);
            W_out_j = [a_j'; W0_j];
        elseif lambda == 0
            X1_j = [ones(1, size(X_j,2)); X_j]; % add intercept
            W_out_j = pinv(X1_j')*Y_j';
            % W_out_j = X1_j' \ Y_j';
        end

        % Residuals
        X_out_j = [ones(1, size(X_j,2)); X_j];
        Y_out_j = W_out_j' * X_out_j;
        Res_j = Y_j - Y_out_j;

        MC(j) = 1 - mean(Res_j.^2);
    end

    MC_by_gamma(i) = min(sum(MC), N);
    fprintf(". Step %d/%d\n", i, L);
end

%% Plot
cmap = turbo(length(gamma_grid)+4);
cmap = cmap(3:end-2,:);

fig = figure(1001);
tiledlayout(1, 1, 'Padding', 'compact')

nexttile
hold on
for l = 1:length(gamma_grid)
    bar(l-1, MC_by_gamma(l), 'FaceColor', cmap(l,:))
end
hold off
set(gca, 'TickDir', 'out');
% xticks(1:2:length(MC_by_gamma))
% xticklabels(arrayfun(@(i) "{10^{"+num2str(i,3)+"}", log10(gamma_grid(1:2:end))))
xt_vals = ceil(log10(gamma_grid_min)*2)/2:0.5:floor(log10(gamma_grid_max)*2)/2;
xt_mapped = (length(gamma_grid)-1) * (xt_vals-log10(gamma_grid_min)) / ...
                (log10(gamma_grid_max)-log10(gamma_grid_min));
xticks(xt_mapped)
xticklabels(arrayfun(@(i) "10^{"+num2str(i,2)+"}", xt_vals))
yticks([1, 10, 20, 30, 40, 50])
yline(N)
yline(1)
ylim([0, N+1])
grid on
xlabel("$\sigma$", 'interpreter', 'latex', 'FontSize', 12)
ylabel("$\widehat{\textnormal{MC}}({\sigma})$", 'interpreter', 'latex','FontSize', 11)

if SAVE_FIGURE
    printpdf(fig, ...
        "./figures/plot_MC_bars_A~"+method_A+".pdf", ...
        [0, 0, 1.1, 1] * 7)
    disp("Figure saved!")
end


% #####
