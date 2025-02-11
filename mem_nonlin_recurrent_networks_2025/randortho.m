function [U] = randortho(N)
%randortho Generate random NxN orthonormal matrix using random normal
%           matrix
% 

[U, ~, ~] = svd(randn(N,N));

end

