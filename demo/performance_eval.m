
%% Performance evaluation
for vv = 1:Realizar
    for zz = 1 :j2
        [~, zMAP]=max(zDist_vv{vv,zz},[],3);
        Z_est_fin1{zz}(:,:,vv) = zMAP;
        labels_Sc_ini{zz}(:,:,vv) = labels_Sc_ini_vv{vv}{zz};
    end
end

%%


%% Shifting of the labels 
% To compute the metrics in the obtained labels considering the
% permutations of the resulting labels

Q_UP = imresize(Z_est_fin1{1},N,'nearest'); % upsampling of the computed labels

lbl = MASK(:);
if lbY == 2
    [Q_UP2,C2_hat] = sw_function(Q_UP,m_ec2,m_ec20,s_ec2,s_ec20,lbl,lbY,Realizar,j2,K,C2);
    %disp('---- C2 - MFA parameter performance-------')
else
    [Q_UP2,~] = sw_function(Q_UP,0,0,0,0,lbl,lbY,Realizar,j2,K,C2);
end


%% analyzing the vectors
%%------- Label  switching--------
tab_to_shufled = [];
lbl = MASK(:);


ec2 = m_ec2;
ec20 = m_ec20;
std_ec2 = s_ec2;
std_ec20 = s_ec20;
for vv1 = 1:Realizar
    A0 = Q_UP(:,:,vv1);
    vec_1 = [mean(3-lbl==(A0(:))),mean(lbl==A0(:))];
    [~,Iu] = max(vec_1);

    if Iu == 1
        Z_estShuffled2(:,:,vv1) = 3-Q_UP(:,:,vv1);
        tab_to_shufled = [tab_to_shufled;vv1];

    else
        Z_estShuffled2(:,:,vv1) = Q_UP(:,:,vv1);

    end
end


if lbY == 2
    % Mean c2 across the realizations
    % C2_real = C2
    set_shuff = unique(tab_to_shufled);
    set_not_shuff = setdiff(1:Realizar,set_shuff)';

    c2_concat   = [ec2(set_not_shuff,:); ec2(set_shuff,[2, 1])];
    c20_concat  = [ec20(set_not_shuff,:); ec20(set_shuff,[2, 1])];

    std_ec2_concat  = [std_ec2(set_not_shuff,:); std_ec2(set_shuff,[2, 1])];
    std_ec20_concat = [std_ec20(set_not_shuff,:); std_ec20(set_shuff,[2, 1])];

    c2_estimado_sw = mean(c2_concat,1);
    c20_f_sw = mean(c20_concat,1);

    c2_f_std_sw = std(c2_concat,[],1);
    c20_f_std_sw = std(c20_concat,[],1);

    std_c2_f_sw = mean(std_ec2_concat,1);
    std_c20_f_sw = mean(std_ec20_concat,1);

    rmse_1_sw = sqrt((c2_estimado_sw-C2).^2 + c2_f_std_sw.^2);
    res_sw = [C2;c2_estimado_sw;c2_f_std_sw;rmse_1_sw;c20_f_sw;c20_f_std_sw];

end



%% Performance evaluation

%--------------------------------------------------------------------------
for vv = 1:Realizar
    Q_aux = Q_UP2(:,:,vv);

    % Evaluation at scale 1
    for kk = 1:K
        A1 = MASK==kk;

        Result = Q_aux==kk;
        [AccF2(vv,kk), SenF2(vv,kk), FmF2(vv,kk), PrF2(vv,kk), MCCF(vv,kk), DiceF2(vv,kk), JacF2(vv,kk), SpeF(vv,kk),errF2(vv,kk)] = EvaluateImageSegmentationScores(A1,Result);


    end
    [oa2(vv), aa2(vv), kappa2(vv), ~, ~] = compute_accuracy(MASK(:),Q_aux(:));

end



disp('Label perfomance mean Dice Coef | Error segmentation')
metrics_ = [mean(DiceF2(:,1),1),mean(DiceF2(:,2),1),mean(errF2(:,1)),mean(errF2(:,2))]
metrics_std = [std(DiceF2(:,1),[],1),std(DiceF2(:,2),[],1),std(errF2(:,1)),std(errF2(:,2))]

%% average to get the error across the regions

errF2(1,:)
mean(DiceF2(1,:),2)
mean(errF2(1,:),2)*100
%--------------------------------------------------------------------------

%%
selected = 1;
figure
imshow(Q_UP2(:,:,selected),[],'Border','tight'),colormap('bone'),axis on
SegError = errF2(selected,:)*100



a_gray = double(MASK);
level = 1;
a_bw = imbinarize(a_gray,level);
Icomp = imcomplement(a_bw);
hold on
[B,L] = bwboundaries(a_bw);


hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'color',[.4 .83 .0], 'LineWidth', 2)
end
hold off
