function  yNZ_k = update_Fourier_Unbiased(K,~,j1,j2,ind0,S_mask_Sc,lx,ind)
%============================================================================
%
%
%
%============================================================================
t = 1; 
I_K = cell(K,1);
for ki = 1: K
    Pv = [];Yv = []; Ni = [];

    for jj=j1:j2
        S = squeeze(S_mask_Sc{jj}(:,:,ki));  
        cte =  sqrt(sum(S(:))); 
        
        %% Fuentes
        Y0 = S.*lx{jj} ;
        Ztilde = sum(Y0(:))/sum(S(:));
        A = Y0 - S.*Ztilde;

        B = fft2(A) ;
        C = fftshift(B);
        Y = C/ cte;

        try
            tp = Y(ind{jj});
            Yv = [Yv tp'];
        catch
            tp=nan(sum(ind{jj}(:)),1);
            Yv = [Yv tp'];
        end
    end
    I.Y(t,:)=Yv;
    I_K{ki} = I;   
end

yNZ_k = cell(K,1);
yZ_k = cell(K,1);
for kj = 1:K
    I_s = I_K{kj};
    [yNZ_k{kj,1},yZ_k{kj,1}] = pyNZ_yZ(I_s.Y,ind0{kj});
end

ny = size(yNZ_k{1},2);      
