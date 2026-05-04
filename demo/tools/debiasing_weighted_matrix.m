function Cgk = debiasing_weighted_matrix(G)
%========================================
% Author: Herwig WENDT
% Function:
%
%========================================
[N1,N2] = size(G);
[rows,cols] = find(G);
ind = sub2ind([N1,N2],rows,cols);
Nt = sum(G(ind));

%% HW 2025/08/25
GG=[G;G;G];
GGG=[GG,GG,GG];
tmp=xcorr2(GGG,G);
Cgk=tmp(N1:N1+2*N1-1,N2:N2+2*N2-1)/Nt;
