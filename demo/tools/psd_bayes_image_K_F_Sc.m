function [LLc,ind0,cov] = psd_bayes_image_K_F_Sc(c2,c20,j1,j2,Nj,ind,Cgn,K)
%==========================================================================
% Debiased Covariance
% Based on the debaised likelihood:
% Paper: The Debiased Spatial Whittle Likelihood Open Access 
%==========================================================================
rc = 3; k = 4;

for kk = 1 : K
    if length(c2) > 1
        c2k = c2(:,kk);c20k = c20(:,kk);
    else
        c20k = c20; c2k = c2;
    end


    Nj1D = Nj(:,1);
    Lc2 = []; Lc20 = [];

    for j=j1:j2

        Cgk = Cgn{j}(:,:,kk);

        % Build full symmetric covariance
        [x,y]  = meshgrid((-(Nj1D(j)):(Nj1D(j)-1)));
        r  = sqrt(x.^2+y.^2);
        rj = floor(Nj1D(j)/k);

        cov.indS{j} = (r<=rc); 
        cov.indL{j} = (r>rc); % Short-covariance

        % Separate contributions of c2 and c20 and store correspondant functions
        cov.var1{j} = log(r/rj);
        rcovL = min(0, (cov.var1{j}).* Cgk ); % *********** here
        rcovL(cov.indS{j}) = 0;
        
        cov.var2{j} = ((log(rc/(rj*2^j))/log(1+rc))*log(1+r)+j*log(2));
        rcovSc2 = (cov.var2{j}).* Cgk ;  % *********** here
        rcovSc2(cov.indL{j})=0;
        
        cov.var3{j} = (1-(log(1+r)/log(4)));
        rcovSc20 = (cov.var3{j}).* Cgk ;  % *********** here
        rcovSc20(cov.indL{j}) = 0;

        % Parameters c2 and c20
        cov.c2{j} = rcovL +rcovSc2;
        cov.c20{j} = rcovSc20;

        % Separate FFT on the two contributions
        tpc2 = real(fftshift(fft2((cov.c2{j}))));
        tpc20 = real(fftshift(fft2((cov.c20{j})))); %FFT

        ic2 = tpc2;
        ic20 = tpc20;

        tpc2 = ic2(2:2:end,2:2:end);
        tpc20 = ic20(2:2:end,2:2:end); % Decimation + vectorization

        % suma en m - low frequencies 
        A = tpc2(ind{j});
        B = tpc20(ind{j});


        Lc2 = [Lc2 A(:)'];
        Lc20 = [Lc20 B(:)'];

        tp(j) = length(A(:));
    end

    ind0{kk} = cumsum([1 tp(j1:j2-1)]);
    LLc{kk}= real([Lc2;Lc20]);



end
