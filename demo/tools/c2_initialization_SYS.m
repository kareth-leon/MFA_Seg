function [Tc2_Est, Tc20_Est] = c2_initialization_SYS(Fi,Gi,iFi,iGi,M,m,ny,K,yNZ_k,LLc,beta2,A)
%==========================================================================
% c2 c20 initalization
%==========================================================================


% Random Values
for kk = 1:K 

    c2m = -0.06; c2M = -0.02;
    c20m = 0.4; c20M = 1.6;

    aa = 0;
    C2i = [c2m + (c2M - c2m ) * rand(M,1) c20m + (c20M-c20m)*rand(M,1)];
    IND = (sum((C2i*LLc{kk}<0),2)>0);

    while sum(IND)>0
        L = length(IND);
        C2new = [c2m+(c2M-c2m)*rand(L,1) c20m+(c20M-c20m)*rand(L,1)];
        C2i(IND) = C2new(IND);IND = (sum((C2i*LLc{kk}<0),2)>0);
        aa = aa + 1;disp(num2str(aa));
    end

    T1i = -C2i(:,1);
    T2i = m{kk}*C2i(:,2) + C2i(:,1);

    T1v(kk) = T1i;
    T2v(kk) = T2i;
end




InnerIter = 50;
for iter = 1: InnerIter
    % if lbY == 2
    for kk = 1 : K
        yNZ = yNZ_k{kk,1}; % Fourier for each region

        % Hidden mean non zero frequency
        T1Fi = T1v(kk)*Fi{kk};
        T2Gi = T2v(kk)*Gi{kk};

        varGG = 1./(1./T1Fi + 1./T2Gi); % ** Eq. (31a): Var Estimation
        muG = varGG.*(1./T1Fi).*yNZ;    % ** Eq. (31a): Var Estimation

        muNZ  = muG + sqrt(varGG/2).*((randn(M,ny) + 1i*randn(M,ny)));

        % v1
        B1 = NormA(yNZ - muNZ, iFi{kk}(:)) + beta2;
        T1v(kk) = (1./gamrnd(A,1./B1));

        % v2
        B2 = NormA(muNZ, iGi{kk}(:)) + beta2;
        T2v(kk) = (1./gamrnd(A,1./B2));

        % Reparametrization of c2 and c^0_2
        Tc2(:,kk,iter) = -T1v(kk);
        Tc20(:,kk,iter)= (T2v(kk) + T1v(kk))/m{kk};
    end

end

inner_burn = round((InnerIter*3)/4);

Tc2_Est = mean(Tc2(:,:,inner_burn:end),3);
Tc20_Est = mean(Tc20(:,:,inner_burn:end),3);