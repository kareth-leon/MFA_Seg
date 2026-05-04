%==========================================================================
% Function to create or load the transition matrices: 
% from parent to current scale (Down)
% child to current scale (Up)
%==========================================================================
if 0
    % Create the transition matrices
    [label_bordeado,pattern_0,n_indjj,TABind,DD_f,DD_t] = Init_Z_Final_Multiscale(labels_Sc,j2,vec_res2,Nj);
    Z_lbl = label_bordeado;
else
    % load the matrices
    if Nj(1,1)==128 && j2==3
        load('ctes_indices_Decimation_matrix.mat','DD_f','DD_t','n_indjj','pattern_0','TABind');
    elseif Nj(1,1)==128 && j2==2
        load('ctes_indices_Decimation_matrix_256_j2_2.mat','DD_f','DD_t','n_indjj','pattern_0','TABind');
    elseif Nj(1,1)==256 && j2==3
        load('ctes_indices_Decimation_matrix_512.mat','DD_f','DD_t','n_indjj','pattern_0','TABind');
    elseif Nj(1,1)==256 && j2==2
        load('ctes_indices_Decimation_matrix_512_j2_2.mat','DD_f','DD_t','n_indjj','pattern_0','TABind');
    end
    Z_lbl{j2,1} = padarray(labels_Sc{j2,1},[1 1],nan,'both');

    for zz = j2-1:-1:1
        Z_lbl{zz,1} = padarray(labels_Sc{zz,1},[vec_res2(zz) vec_res2(zz)],nan,'both');
    end
end
