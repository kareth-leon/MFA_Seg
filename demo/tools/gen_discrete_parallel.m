function Y = gen_discrete_parallel(X,p,N,M)

p = p./(ones(size(p,1),1)*sum(p,1));
cum_p = [zeros(1,size(p,2)); cumsum(p,1)];

NM = N*M;

z = rand(1,NM);
Y(NM)=0;
for i=1:NM
    ind = find(z(i)<cum_p(:,i), 1 );
    Y(i) = X(ind-1);
end
 
Y = reshape(Y,N,M);