function [mm,Fi,Gi,iFi,iGi] = psd_bayes_image_K(j1,j2,cov,ind,Cgn,K,ind0)
% Ligther version 
%========================================
% Function: debiased covariance
%
%========================================
rc = 3; 


Fi = cell(K,1);
iFi = cell(K,1);
Gi = cell(K,1);
iGi = cell(K,1);
mm = cell(K,1);


Av  = cov.var1;
indS= cov.indS;
Av2 = cov.var2;

indL= cov.indL;
Av3 = cov.var3;

for kk = 1 : K
    Lc2 = [];Lc20 = [];    
    for j=j1:j2

        Cgk = Cgn{j}(:,:,kk);
        
        % Separate contributions of c2 and c20 and store correspondant functions
        rcovL = min(0, (Av{j}).* Cgk ); % *********** here
        rcovL(indS{j}) = 0;
        
        rcovSc2 = (Av2{j}).* Cgk ;  % *********** here
        rcovSc2(indL{j})=0;
        
        rcovSc20 = (Av3{j}).* Cgk ;  % *********** here
        rcovSc20(indL{j}) = 0;

        tpc2 = real(fftshift(fft2((rcovL +rcovSc2))));
        tpc20 = real(fftshift(fft2((rcovSc20)))); %FFT

        ic2 = tpc2;
        ic20 = tpc20;

        tpc2 = ic2(2:2:end,2:2:end);
        tpc20 = ic20(2:2:end,2:2:end); % Decimation + vectorization


        A = tpc2(ind{j});
        B = tpc20(ind{j});

        Lc2 = [Lc2 A(:)'];
        Lc20 = [Lc20 B(:)'];

        tp(j) = length(A(:));
        
    end



    %% pLLc_prior_K

    auxLLc =real([Lc2;Lc20]);

    indM = find(-auxLLc(1,:)<0);
    m  = min(-auxLLc(2,indM)./(-auxLLc(1,indM)))-10^(-3.5);
    mm{kk}  = m;

    % Non-zero frequency
    LLcNZ = auxLLc;
    LLcNZ(:,ind0{kk}) = [];

    fi  = -LLcNZ(1,:);
    gi = LLcNZ(2,:);

    Fi{kk}  = (m * fi + gi )/m;
    Gi{kk}  = gi / m;

    iFi{kk} = 1./Fi{kk};    
    iGi{kk} = 1./Gi{kk};

end
