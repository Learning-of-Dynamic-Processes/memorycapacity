function [MCs, Q, h] = mc_arnoldi(A, C, l)
%MC_ARNOLDI Compute the LESN memorie using Arnoldi iterations
%

N = size(A,1);

if nargin == 2
    l = N;
else
    if l > N
        l = N;
    end
end

Q = zeros(N,l);
h = zeros(N,l+1);
q = C / norm(C); % normalize
Q(:,1) = q;
for k = 1:(l-1)
    v = A * q;
    for j = 1:k
        h(j,k) = dot(conj(Q(:,j)), v);
        v = v - h(j,k)*Q(:,j);
    end
    h(k+1,k) = norm(v);
    if h(k+1,k) > 10^-13
        q = v / h(k+1,k);
        Q(:,k+1) = q; 
    else
        break
    end
end

% Compute MCs
MCs = diag(Q' / (Q * Q' + 10^-14*eye(N)) * Q);

end

