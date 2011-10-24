
J = 5^2;
L = 4;
M = 3;
T = 2;
lambda = 1;

c = randn(M, 1);
psi = randn(L, L, M);
X0 = randn(L, 1);
X1 = randn(L, 1);
V = randn(J,L);


Jsz = sqrt(J);   
border = 1;
mask = zeros(Jsz,Jsz);
valid_range = border+1:Jsz-border;
mask(valid_range,valid_range) = ones(length(valid_range));
mask = reshape(mask, J, 1);


tic
checkgrad('objfun_c', c(:), 1e-4, psi, X0, X1, V, lambda, mask)
toc

