function [W_out, Res_Y, Y_out, Y, X] = ESN(Yt, zt, A, C, fun, lambda, discard)
%ESN Estimate a linear ESN model
%   
T = size(Yt, 2);
K = size(A, 1);

% State
X = zeros(K, T);
X(:,1) = C*zt(:,1);
for t = 2:T
    X(:,t) = fun(A*X(:,t-1) + C*zt(:,t));
end
X = X(:,(discard+1):end);
Y = Yt(:,(discard+1):end);

% Regression
if lambda > 0
    W0 = (X*X' + T*lambda*eye(K)) \ (X*Y');
    a = mean(Y,2) - W0' * mean(X,2);
    W_out = [a'; W0];
elseif lambda == 0
    X1 = [ones(1, T-discard); X]; % add intercept
    W_out = pinv(X1')*Y';
end

% Rss = @(W) sum((Y - W * X).^2);
% W_out = fminunc(Rss, [1, zeros(1, K)])';

X = [ones(1, T-discard); X];

Y_out = W_out' * X;
Res_Y = Y - Y_out;

end