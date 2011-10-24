function [f,g] = objfun_c(x0,psi,X0,X1,V,lambda,mask);

L = size(X0,1);
Lsz = sqrt(L);

c = x0;

M = size(c, 1);

A = zeros(L);

for i = 1:M
    A = A + psi(:,:,i) * c(i);
end

[U,D] = eig(A);
W = inv(U).';
D = diag(D);

ExpA = real(expm(A));


%E = X1 - ExpA*X0;
%E = V*(X1 - ExpA*X0);

E = mask .* (V*(X1 - ExpA*X0));

f = 0.5*sum(E(:).^2) + lambda*sum(abs(c(:)));

%% populate F
F = zeros(L,L);

ExpD = exp(D);
for i = 1:L
    for j = 1:L
        if D(i) == D(j)
            F(i,j) = ExpD(i);
        else
            F(i,j) = (ExpD(j) - ExpD(i)) / (D(j) - D(i));
        end
    end
end


psi_2d = reshape(psi, L^2, M);

VtV = V'*diag(mask)*V;
dc_exp = -W * (F .* (U.' * ( VtV*X1*X0' - VtV*ExpA*X0*X0' ) * W)) * U.';

if 1
    dc = dc_exp(:)'*psi_2d;
    dc = dc + lambda*sign(c)';
else
    dc = zeros(M,1);
    for i = 1:M
        dc(i) = sum(sum( dc_exp .* psi(:,:,i) ));
    end
    dc = dc + lambda*sign(c);
end

g = real(dc(:));

if ~isfinite(f) || sum(~isfinite(g)) > 0
    g = zeros(size(g));
    f = sqrt(realmax);
end

