function [f,df] = objfun_psi(x0,c,X0,X1,V,mask);

L = size(X0,1);
M = size(c,1);

psi = reshape(x0, L, L, M);

A = zeros(L);
for i = 1:M
    A = A + psi(:,:,i) * c(i);
end

[U,D] = eig(A);
W = inv(U).';

D = diag(D);


ExpA = U * diag(exp(D)) * W.';

%E = X1 - ExpA*X0;
E = mask .* (V*(X1 - ExpA*X0));


f = 0.5*sum(E(:).^2);
f = real(f);


if nargout > 1

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


    VtV = V'*diag(mask)*V;
    dpsi_exp = -W * (F .* (U.' * ( VtV*X1*X0' - VtV*ExpA*X0*X0' ) * W)) * U.';
    dpsi = dpsi_exp(:)*c';

    df = real(dpsi(:));

end


