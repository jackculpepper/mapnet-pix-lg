
load_X

psi = randn(Lmax,Lmax,N);
for n = 1:N
    psi(:,:,n) = sqrt(Lmin)*psi(:,:,n)/sqrt(sum(sum(psi(:,:,n).^2)));      
end

update = 1;

