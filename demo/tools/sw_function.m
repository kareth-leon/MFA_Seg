function [Z_estShuffled2,res_sw] = sw_function(X,ec2,ec20,std_ec2,std_ec20,MASK,lbY,Realizar,j2,K,C2)


%% Label  switching
tab_to_shufled = [];
lbl = MASK(:);
if K ==2

    for vv1 = 1:Realizar
        A0 = X(:,:,vv1);
        vec_1 = [mean(3-lbl==(A0(:))),mean(lbl==A0(:))];
        [~,Iu] = max(vec_1);

        if Iu == 1
            Z_estShuffled2(:,:,vv1) = 3-X(:,:,vv1);
            tab_to_shufled = [tab_to_shufled;vv1];

        else
            Z_estShuffled2(:,:,vv1) = X(:,:,vv1);

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

else

    vec_permut = circshift(perms(1:K),1);
    tab_to_shufled_data = [];
    factor_shift = [];

    aux_sc = 1;

    for vv1 = 1:Realizar

        A0 = Q_UP(:,:,vv1);
        % original state
        for kk = 1:K
            A00(:,:,kk) = A0 == kk;
        end

        vec_vec_1 = [];
        vec_vec_2 = [];

        sz = size(vec_permut,1);

        % ------------------------------------------------------%
        % ------ prueban todas las posibles combinaciones-------%
        % ------------------------------------------------------%
        for vvo = 1: sz
            rt = vec_permut(vvo,:);
            [~,perm_A00{vvo}] = max(A00(:,:,rt),[],3);
            vec_vec_1 = [vec_vec_1,mean(lbl(:)==perm_A00{vvo}(:))];
            if 0
                for kk = 1:K
                    A1 = label_mtx01{aux_sc}==kk;
                    B1 = perm_A00{vvo}==kk;
                    [Acc(kk), Sen(kk), Fm(kk), Pr(kk), MCC(kk), Dice(kk), Jac(kk), Spe(kk)] = EvaluateImageSegmentationScores(A1,B1);
                end
                vec_vec_2 = [vec_vec_2,mean(Spe(:))];
            else
                [oa, aa, kappa, class_acc, confusion] = compute_accuracy(lbl(:),perm_A00{vvo}(:));
                vec_vec_2 = [vec_vec_2,kappa];
            end
        end

        clear A00
        % vector con el promedio mas alto de que sean iguales al
        % groundtruth
        [a,Iu] = max(vec_vec_2);

        rt1 = vec_permut(Iu,:);
        factor_shift = [factor_shift;rt1];

        if exist("A001",'var'); clear A001;end

        Akk = Q_UP(:,:,vv1);
        for kk = 1:K
            % 1. se busca donde están las posiciones originales
            A001(:,:,kk) = Akk == kk;
        end
        [~,last11] = max(A001(:,:,rt1),[],3);
        Q_UP2(:,:,vv1) = last11;

        if Iu ~= sz +1 % no shufling
            tab_to_shufled_data = [tab_to_shufled_data;vv1];
        end

        clear A001;


    end
    if lbY == 2
        c2_concat = [];
        c20_concat = [];
        std_ec2_concat = [];
        std_ec20_concat = [];
        for ii = 1:size(factor_shift,1)
            c2_concat = [c2_concat;ec2(ii,factor_shift(ii,:))];
            c20_concat= [c20_concat;ec20(ii,factor_shift(ii,:))];

            std_ec2_concat= [std_ec2_concat;std_ec2(ii,factor_shift(ii,:))];
            std_ec20_concat= [std_ec20_concat;std_ec2(ii,factor_shift(ii,:))];
        end

        c2_estimado_sw = mean(c2_concat,1);
        c20_f_sw = mean(c20_concat,1);

        c2_f_std_sw = std(c2_concat,[],1);
        c20_f_std_sw = std(c20_concat,[],1);

        std_c2_f_sw = mean(std_ec2_concat,1);
        std_c20_f_sw = mean(std_ec20_concat,1);

        rmse_1_sw = sqrt((c2_estimado_sw-C2).^2 + c2_f_std_sw.^2);
        C2_hat = [C2;c2_estimado_sw;c2_f_std_sw;rmse_1_sw;c20_f_sw;c20_f_std_sw];

    else
        C2_hat = [];
    end

    %
    % else
    %     Z_estShuffled2 = X ;
    %     res_sw =  [];

end