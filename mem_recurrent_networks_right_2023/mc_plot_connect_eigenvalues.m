%
%
%

clear
rng('default')

% Save plots?
SAVE_PLOT_FLAG = false;

%% Parameters

N = 100;

% Method to sample the connectivity matrix, spectral radius
method_A = 'sparse_normal_cond';

switch method_A
    case 'normal'
        A = randn(N) / sqrt(N);
    case 'uniform'
        A = (rand(N)-.5) * 2 / sqrt(N/3);
    case 'sparse_normal'
        A = full(sprandn(N,N,0.1)) / sqrt(N*0.1);
    case 'sparse_normal_cond'
        A = full(sprandn(N,N,0.1,0.7));
    case 'ortho'
        A = randortho(N);
    case 'cyclic'
        A_ = eye(N);
        A = A_(:,[N, 1:N-1]);
    case 'delay'
        A = toeplitz([0, 1, zeros(1,N-2)], zeros(1,N));
    otherwise
        error("No appropriate sampling rule for A chosen.")
end

%% Plot

fig = plotConnectEigs(A);
% hold on
% fimplicit(@(x,y) x.^2 + y.^2 - rho);
% hold off
legend off 
% fig = plotConnectSvals(A);

% Save figure
if SAVE_PLOT_FLAG
printpdf(fig, ...
    join(["figures/", "plot_A_eigs", ...
            "_N=", string(N), ...
            "_A=", string(method_A), ".pdf"],""), ...
    [0, 0, 8, 6])
end

% #####