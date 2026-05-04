%% Initial Parameters
%--------------------------------------------------------------------------
zDist = cell(j2,1);
for zz = 1:j2
    Njjj = Njj(zz,:);
    zDist{zz,1} = zeros([Njj(zz,:),K]);
    NjjjNUEVA{zz} = size(pattern_0{zz});
    mtx_Z_iterc{zz,1}(:,:,1) = Z_lbl{zz,1};

end

%% Parameters labels
% Constant parameters  in the label estimation
pmlbl.DD_t = DD_t;
pmlbl.DD_f = DD_f;
pmlbl.NjjjNUEVA = NjjjNUEVA;
pmlbl.TABind = TABind;
pmlbl.n_indjj=n_indjj;
pmlbl.j2=j2;
pmlbl.K=K;
pmlbl.beta_current = beta_current;
beta_current_init = beta_current;

if plot1 == 1
fig = figure;
colormap(flipud(gray));
set(fig, 'Position', [100, 100, 1000, 500]); % Define figure size early on
end

iter2=2;
iter3=2;
for iter = 1 :NM
    % iter
    %------------------------------------------------------
    %                           MFA PARAMETERS
    %------------------------------------------------------
    if lbY == 2
        % ==============================================================
        % This part is based on Bayesian estimation for multifractal analysis
        %
        % Paper: Wendt, Herwig, et al. "Multifractal analysis of multivariate 
        % images using gamma Markov random field priors." SIAM Journal on 
        % Imaging Sciences 11.2 (2018): 1294-1316.
        %
        % Toolbox: https://www.irit.fr/~Herwig.Wendt/software.html#bayesc2
        % ==============================================================

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
            Tc2(:,kk,iter+1) = -T1v(kk);
            Tc20(:,kk,iter+1)= (T2v(kk) + T1v(kk))/m{kk};
        end
    end

    %------------------------------------------------------
    %                         LIKELIHOOD
    %------------------------------------------------------
    % likelihood estimation
    for zz = j1  : j2
        Njjj = Nj(zz,:);
        for kk = 1: K
            if lbY == 2
                c02 = Tc20(:,kk,iter+1);c2= Tc2(:,kk,iter+1);
            else
                c02=c_20_est2(kk);c2=C22(kk);
            end
            part_1 = sqrt(c02 + c2*log(2^(zz)));

            % Log-Likelikhood
            loglike{zz,1}(:,:,kk) = reshape(log(normpdf(xx{zz},0,(part_1))),Njjj);

            loglike_border{zz,1}(:,:,kk) = padarray(loglike{zz,1}(:,:,kk),[vec_res2(zz) vec_res2(zz)],nan,'both');
        end
    end


    %------------------------------------------------------
    %                         LABELS
    %------------------------------------------------------
    Z_lbl = updateChessboard_MATLAB_SYS(Z_lbl,loglike_border,pmlbl);



    %
    % %------------------------------------------------------
    % %                         show
    % %------------------------------------------------------
    if plot1 == 1

        subplot(3,j2,1)
        imshow(Z_lbl{1},[]),title('j_1')

        subplot(3,j2,2)
        imshow(Z_lbl{2},[]),title('j_2')

        if j2 == 3
            subplot(3,j2,3)
            imshow(Z_lbl{3},[]),title('j_3')
        end
        if lbY == 2
            subplot(3,j2,[4,6])
            plot(squeeze(Tc2(:,1,:)))
            title(['Theta 1: ' num2str(Tc2(:,1,iter))]),grid on

            subplot(3,j2,[7,9])
            plot(squeeze(Tc2(:,2,:))),xlabel('Iteration'),
            title(['Theta 2: ' num2str(Tc2(:,2,iter))]),grid on
            if K == 3
                subplot(3,j2,6)
                plot(squeeze(Tc2(:,3,:)))%,xlabel('Iter'),
                title(num2str(Tc2(:,3,iter))),grid on
            end
        end
        sgtitle(['Iteration: ' num2str(iter)])
        drawnow
    end

    %----------------------------------------------------------------------
    %************* CHECK ALL CLASSES WITH VALUES ****** & SAVE EACH ITER
    %----------------------------------------------------------------------
    for zz = 1:j2
        zzv = 0;
        for kk = 1: K
            idx = Z_lbl{zz,1} == kk;
            if~(sum(idx(:)))
                zzv = 1 ; %[zzv,zz];
                disp('--------------------------------------------------------')
                disp(['collapseee at iteration: ' num2str(iter) ', zz:' num2str(kk) ', kk:' num2str(zz)])
                disp('--------------------------------------------------------')
            end
        end
        if zzv == 1
            Z_lbl{zz,1} = mtx_Z_iterc{zz,1}(:,:,iter) ; 
        end
        mtx_Z_iterc{zz,1}(:,:,iter+1) = Z_lbl{zz,1};

        %------------------------------------------------------
        % Optimizando este paso para no hacerlo dos veces
        for kk = 1:K
            Z_est_cut{zz,1}(:,:,kk) = Z_lbl{zz,1}(vec_res2(zz)+1:end-vec_res2(zz),vec_res2(zz)+1:end-vec_res2(zz))==kk;
        end
    end

    %------------------------------------------------------
    %   UPDATE MAP
    %------------------------------------------------------
    if iter >= burnin
        for zz = 1:j2
            zDist{zz,1} = zDist{zz,1} + Z_est_cut{zz,1};
            clear aux_Z
        end
    end

    %------------------------------------------------------
    %                         UPDATE VARIABLES
    %------------------------------------------------------
    if lbY == 2
        % ------------------------------------------------------------
        % #2: UPDATE: Fourier coefficients given the new MASK
        yNZ_k = update_Fourier_Unbiased(K,Nj,j1,j2,ind0,Z_est_cut,lx,indd);


        % -------------------------------------------------------------
        % #3: UPDATE: Cgn given the new MASK
        for j = j1: j2
            for kk = 1 :K
                Cgn_Sc{j}(:,:,kk) = debiasing_weighted_matrix(squeeze(Z_est_cut{j}(:,:,kk)));

            end
        end
        % -------------------------------------------------------------
        % #4: UPDATE: LLc and then Fi Gi and m
        [m,Fi,Gi,iFi,iGi] = psd_bayes_image_K(j1,j2,cov,indd,Cgn_Sc,K,ind0);
        % -------------------------------------------------------------
    end

    %------------------------------------------------------
    %                         UPDATE BETA
    %------------------------------------------------------
    if SB_Y == 1
        if (iter <= burnin)
        % #5: UPDATE: beta parameter
            for more_iter = 1:pmlbl.it2

                a_nt = 10*(prod(Njj').^-1)*(iter2-1)^(-3/4);

                betaNew = updateGranularity(Z_lbl,Z_est_cut,pmlbl,loglike_border,vec_res2,Njj,iter2);

                pmlbl.beta_current = betaNew;

                beta_iter(:,:,iter3) = betaNew;

                iter2 = iter2+1;

            end

            iter3 = iter3 + 1; % for the beta tensor across iterations!

        end
    end
end


% Save
if lbY == 2
    output.Tc2 = Tc2;
    output.Tc20 = Tc20;
end


output.zDist = zDist;
% output.beta_pott = beta_pott;
% output.time1 = time1;
vv = 1;
if exist('beta_iter','var')==0
    beta_iter = [];
end

output.beta_iter = beta_iter;



output.mtx_Z_iterc = mtx_Z_iterc;


if lbY == 2
    c2_Realizar{vv} = output.Tc2;
    c20_Realizar{vv} = output.Tc20;

    % MMSE ~ average
    for kk=1:K
        m_ec2(vv,kk) = mean(output.Tc2(:,kk,burnin:NM),3);
        m_ec20(vv,kk) = mean(output.Tc20(:,kk,burnin:NM),3);

        s_ec2(vv,kk) = std(output.Tc2(:,kk,burnin:NM),[],3);
        s_ec20(vv,kk) = std(output.Tc20(:,kk,burnin:NM),[],3);
    end
end

for zz = 1 :j2
    zDist_vv{vv,zz} = output.zDist{zz};

end
mtx_Z_iterc2{vv} = output.mtx_Z_iterc;
labels_Sc_ini_vv{vv} = output.labels_Sc_ini;

if SB_Y == 1
    beta_real{vv} = output.beta_iter;
end
