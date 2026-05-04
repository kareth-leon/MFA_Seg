function [lx,vec_feat,Ni,Nj,I_K,Lx_Z_k,I_ori,S_mask_Sc,Cgn_Sc,labels_Sc] = pDWT_2DNZ_SYS(data,j1,j2,Nwt,gamint,K,label_mtx01)
%==========================================================================
% Log-leaders estimation and labels initialization using the logleaders
%==========================================================================

%% Wavelet leaders
j2_param = j2;

eta = 0.25;
t = 1;
tmp1=data(:,:,t);

if sum(~isfinite(tmp1(:)))
    for j = j1: j2
        leaders(j).matrix=inf(size(data(1:2^j:end,1:2^j:end,t)));
    end
else
    % Wavelet transform
    [~, leaders, ~] = DxLx2d_bayes_dwt(data(:,:,t),Nwt,j2,gamint);
    if length(leaders)<j2 && ~sum(isfinite(leaders(end).allmatrix(:)))
        for j=length(leaders)+1:j2
            leaders(j).allmatrix = leaders(j-1).allmatrix(1:2:end,1:2:end);
        end
    end
end
if length(leaders)<j2
    fprintf('\n\n ERROR:\n There are less than j2 scales available.\n Reduce j2 and restart.\n\n');
    error(' ');
end

for j=j1:j2
    lx{j}   = log(leaders(j).matrix);
    Lx{j}{t} = lx{j};
    VV(t,j) = var(lx{j}(:));
    mV(t,j) = mean(lx{j}(:));
    Nj(j,:) = size(lx{j});
end

[~,ind] = fft_freqNZ(j1,j2,Nj,eta);


%% Labels initialization

if 1
    %% cutting the images
    if j2 == 5
        j21 = j2-1;
    else
        j21 = j2;
    end
    OL = 0.75;

    if j2 == 2;patch_size=[16,8];end
    if j2 == 3;patch_size=[16,8,4];end
    if j2 == 4;patch_size=[32,16,8,4];end
    if j2 == 5;patch_size=[32,16,8,4,2];end

    % build set of patches
    for zz = j1:j21
        sw = patch_size(zz)*[1 1];
        [log_Data{zz},pdata{zz},u{zz},b0{zz}] = cuteImage(lx{zz},sw,OL*(sw));
    end

    M = length(log_Data{1});

    for t=1:M
        for zz=j1:j21
            lx0{zz}   = log_Data{zz}(:,:,t);
            %Lx{zz}{t} = lx0{zz};
            VVini(t,zz) = var(lx0{zz}(:));
            VV2ini(t,zz) = std(lx0{zz}(:));
            mVini(t,zz) = mean(lx0{zz}(:));
            mV2ini(t,zz) = median(lx0{zz}(:));
            Njini(zz,:) = size(lx0{zz});
            min_(t,zz) = min(lx0{zz}(:));
            max_(t,zz) = max(lx0{zz}(:));
        end
    end


    [ec2,ec20] = regr_cp(VVini,j1,j21,Njini(:,2)',0,[sqrt(length(VVini)),sqrt(length(VVini))] );

    [ec1,ec10] = regr_cp(mVini,j1,j21,Njini(:,2)',0,[sqrt(length(VVini)),sqrt(length(VVini))] );


    sigma = 2; %[0.1 0.5 1];
    vec_sig =[];
    vec_sig2 =[];

    % C2, C20
    for fil = 1:length(sigma)
        aux = imgaussfilt(ec2,sigma(fil));
        aux2 = imgaussfilt(ec20,sigma(fil));

        aux1 = imgaussfilt(ec1,sigma(fil));
        aux3 = imgaussfilt(ec10,sigma(fil));

        c22 = [aux(:),aux2(:)];
        vec_sig = [vec_sig,c22];

        vec_sig2 = [vec_sig2,c22,aux1(:),aux3(:)];
    end

    vec_feat{1} = ec2(:);
    vec_feat{3} = ec20;

    vec_feat{7} = aux(:);
    vec_feat{2} = aux2(:);
    name_2{2} = {'3:[fil c20]'};

    vec_feat{6} = [min_,max_];
    name_2{6} = {'6:[min_,max_]'};

    vec_feat{8} = [ec2(:),ec20(:)];
    vec_feat{9} = c22;
    vec_feat{10} = [ec20(:),aux2(:)];


    g = 2;
    idx1 = kmeans(vec_feat{g},K,'Replicates',10);%,'start','cluster');%);%'

    idx_R{g} = reshape(idx1,[sqrt(length(VVini)),sqrt(length(VVini))]);


    label_mtx{j2} =  imresize(idx_R{g},Nj(j2,:),'nearest');
    for jj=j2-1:-1:1
        aux_lbl = label_mtx{jj+1};
        label_mtx{jj} = couting_spatial_replicate(aux_lbl);
    end
else
    for zz = 1:j2
        Nnj = Nj(zz,:)/2;
        label_mtx{zz,1} = randi([1,K],[Nj(zz,:)]);
    end


    vec_feat = [];
end


%% 
for j=j1:j2
    for ki = 1:K
        aux = label_mtx{j};
        S_mask(:,:,ki) = aux ==ki;
        S_mask_Sc{j,1}(:,:,ki) = S_mask(:,:,ki);
        Cgn_Sc{j}(:,:,ki) = debiasing_weighted_matrix(S_mask(:,:,ki));
    end
    clear S_mask
    labels_Sc{j,1} = aux;
end



t = 1;
%% Function fourier

I_K = cell(K,1);
for ki = 1: K
    Pv = [];Yv = []; Ni = [];

    for j=j1:j2_param

        S = S_mask_Sc{j}(:,:,ki);
        cte =  sqrt(sum(S(:))); % fijo

        %% Fuentes
        Y0 = S.*lx{j};
        Ztilde = sum(Y0(:))/sum(S(:));
        Lx_Z_k{j}(:,:,ki) = Y0 - S.*Ztilde;

        B = fft2(Lx_Z_k{j}(:,:,ki)) ;
        C = fftshift(B);
        Y = C/ cte;

        P = abs(Y).^2;
        try
            tp = P(ind{j});
            Pv = [Pv tp'];
        catch
            tp=nan(sum(ind{j}(:)),1);
            Pv = [Pv tp'];
        end
        try
            tp = Y(ind{j});
            Yv = [Yv tp'];
        catch
            tp=nan(sum(ind{j}(:)),1);
            Yv = [Yv tp'];
            disp('entró')
        end
    end

    I.Y(t,:)=Yv;
    I.P(t,:)=Pv;

    I_K{ki} = I;

end



%% Fourier transform

Pv = [];Yv = []; Ni = [];

for j=j1:j2_param
    C = fftshift(fft2(lx{j}));
    Y = C / sqrt(prod(Nj(j,:)));

    P = abs(Y).^2;
    try
        tp = P(ind{j});
        Pv = [Pv tp'];
    catch
        tp=nan(sum(ind{j}(:)),1);
        Pv = [Pv tp'];
    end
    try
        tp = Y(ind{j});
        Yv = [Yv tp'];
    catch
        tp=nan(sum(ind{j}(:)),1);
        Yv = [Yv tp'];
    end
end
if t==1
    Ni = ones(size(Yv));
end
I_ori.Y(t,:)=Yv;
I_ori.P(t,:)=Pv;



end



