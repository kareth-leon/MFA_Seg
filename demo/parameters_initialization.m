%% Default parameters in Paper simulations
%--------------------------------------------------------------------------
% Experiment's parameters
NM      = 300; % Number of iterations
Realizar= 1; % Number of realizations
lbY     = 2; % Joint estimation: labels and params? Y:2, N:1
K       = 2; % Number of multifractal values in the image
SB_Y    = 1; % Y:1,N:0
save_if = 1; % save all
plot1   = 1;
number_run = 1; % number of run with all the same

Nwt     = 1;    % number of vanishing moments of wavelet (Daubechies)
j1      = 1;    % scaling range lower cutoff
j2      = 3;    % multifractal analysis parameters

%--------------------------------------------------------------------------
gamint  = 1;    % fractional integration parameter for wavelet leaders

% Initial beta vector for K = 2
beta_1   = [4,0.5,4,0.5,0.5,0.5]; % beta values
pmlbl.it = 2;   % Iterations to update beta
pmlbl.it2= 1;   % Iterations to update beta

%% Initialize multifractal parameters
if ~exist('C2','var'); C2 = zeros(1,K); end
c_20_est = zeros(1,K);


%% params
[mm,nn,tt] = size(data);
N = [mm,nn];
%--------------------------------------------------------------------------
bi = 10; % Percentage of samples for burnin
burnin = round(NM*(bi/100));
if burnin == 0; burnin = 1;end
%--------------------------------------------------------------------------
% Compute Fourier coefficients for the Whittle approximation
[Njj,indd,~,~,~] = pDWT_2DNZ(randn(N(1),N(2)),j1,j2,Nwt,0.25);
%--------------------------------------------------------------------------

%% Create mask for the different scales
for zz = 1:j2
    label_mtx0 = imresize(MASK,Njj(zz,:),'nearest');
    label_mtx01{zz,1} = label_mtx0;
    groundtruth{zz,1} = label_mtx0(:);
end

%% params for the beta estimation
%--------------------------------------------------------------------------
beta_base = zeros(3,3);
[bn,bm] = size(beta_1);
bi_vec  = [];
bi_r    = 1;

beta_current = beta_base;
beta_current(1,:) = beta_1(bi_r,1:3);% spatial
beta_current(3,:) = beta_1(bi_r,4:6); % scales
beta_iter(:,:,1) = beta_current; % save the computed betas

beta_1 = [3,1,3,0.5,0.5,0.5];
pmlbl.it = 5;   % Iteraciones dentro del update de beta
pmlbl.it2= 2;   % Iteraciones fuera del update de beta

%% Compute the logleaders coefficients
[lx,~,~,Nj,I_K_ini,~,~,~,Cgn_Sc,labels_Sc0] = pDWT_2DNZ_SYS(X,1,j2,Nwt,gamint,K);

%------------------------------------------------------------------
output.lx = lx;
labels_Sc = labels_Sc0;

for zz = j1:j2
    xx{zz} = lx{zz}-mean(lx{zz}(:));
    labels_Sc_ini{zz} = labels_Sc0{zz};

end
output.labels_Sc_ini = labels_Sc_ini;
Cgn_Sc_ini = Cgn_Sc;

%------------------------------------------------------------------

if lbY == 1
    a_pixel = labels_Sc0{j1}(12,12); % 1
    b_pixel = labels_Sc0{j1}(64,64); % 2
    if a_pixel~=1 && b_pixel ~= 2
        c_20_est2 = c_20_est(end:-1:1);C22 = C2(end:-1:1);
    else
        c_20_est2 = c_20_est;C22 = C2;
    end
else
    c_20_est2 = c_20_est;C22 = C2;
end

%--------------------------------------------------------------------------
% COMPUTE PSD VIA FFT OF THE COVARIANCE
if lbY == 2

    [LLc,ind0,cov] = psd_bayes_image_K_F_Sc(-0.04,0.6,j1,j2,Njj,indd,Cgn_Sc_ini,K);
    
    [m,Fi,Gi,iFi,iGi,pL] = pLLc_prior_K(LLc,j1,j2,Njj,ind0,0,K);

    % Fourier estimations
    [ny,yNZ_k, M ] = fourier_est(I_K_ini,ind0,K);
    M2=M;

    yNZ_k2 = yNZ_k;

    alpha = 1e-3;
    beta2 = 1e-3;
    A = ny+alpha;

    %--------------------------------------------------
    % Initialization of the parameters MFA
    [Tc2_ini, Tc20_ini] = c2_initialization_SYS(Fi,Gi,iFi,iGi,M,m,ny,K,yNZ_k,LLc,beta2,A);

    disp('Initial C2')
    disp(num2str(Tc2_ini))

    disp('Initial C20')
    disp(num2str(Tc20_ini))

    %--------------------------------------------------
    for kk = 1: K
        Tc2(:,kk,1)   = Tc2_ini(:,kk);
        Tc20(:,kk,1)  = Tc20_ini(:,kk);

        %--------------------------------------------------
        % Vector initialization
        % for potts sampling
        T1v(kk) = -Tc2(:,kk);
        T2v(kk) = m{kk}*Tc20(:,kk) + Tc2(:,kk);

    end
end
%--------------------------------------------------------------------------
% Labels init value.
lbl_temp = labels_Sc;
lbl_temp2 = labels_Sc;

% Scale relationship
vec_res = [1,2,4,8,16];
vec_res2 = vec_res(j2:-1:1); 