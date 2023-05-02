%
%
%

clear
rng('default')
 
% Save plots?
SAVE_PLOT_FLAG = false;

%% Parameters

N = 50;

% Method to sample the connectivity matrix, spectral radius
method_A = {'normal', 'sparse_normal', 'uniform', 'ortho', 'cyclic'};
rho = 0.9;

% Method to sample the input mask
method_C = 'normal';

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

store = cell(length(method_A),1);
for m = 1:length(method_A)
    methodA = method_A{m};
    
    switch methodA
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
    
%     X = zeros(N, max(1000, 5*N));
%     Ac = eye(N)*C;
%     for i = 1:size(X,2)
%         X(:,i) = Ac;
%         Ac = A * Ac;
%     end
%     Gamma = X * X';

    Gamma = C * C';
    P = A;

    while norm(P) > eps
        Gamma = P * Gamma * P' + Gamma;
        P = P * P;
    end
    
%     for i = 1:20
%         Gamma = P * Gamma * P' + Gamma;
%         P = P * P;
%     end
    
    spec_Gamma = sort(abs(eig(Gamma)), 'desc');
    
    store{m}.methodA = methodA;
    store{m}.specGamma = spec_Gamma;
end

%% Plot

fig = figure();
semilogy(1:N, store{1}.specGamma)
hold on
for m = 2:length(store)
    semilogy(1:N, store{m}.specGamma)
end
yline(eps, 'k-', "eps", 'LabelHorizontalAlignment', 'right')
grid
title("Eigenvalues of $G_{\textnormal{x}}$", 'interpreter', 'latex')
legend(cellfun(@(x) strrep(string(x.methodA), "_", " "), store), ...
        'Location', 'best')

% Save figure
if SAVE_PLOT_FLAG
printpdf(fig, ...
    join(["figures/", "plot_gamma_eigs", ...
            "_N=", string(N), ...
            ".pdf"],""), ...
    [0, 0, 12, 10])
end

% #####