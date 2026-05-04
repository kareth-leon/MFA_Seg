function [z_bordeado,pattern_0,n_indjj,TABind,DD_f,DD_t] = Init_Z_Final_Multiscale(labels_Sc,j2,vec_res2,Nj)
%==========================================================================
% create decimation matrices when N =256, j2=3
% chess board multiscale
%
%==========================================================================


v = ones(size(labels_Sc{j2,1})); v(2:2:end)=0;
to_2 = toeplitz(v);

pattern_0{j2,1} = padarray(to_2(1:Nj(j2,1),1:Nj(j2,1)),[1 1],nan,'both');
z_bordeado{j2,1} = padarray(labels_Sc{j2,1},[1 1],nan,'both');
for zz = j2-1:-1:1
    pattern_0{zz,1} = repelem(pattern_0{j2,1},vec_res2(zz),vec_res2(zz));
end

%%
for zz = j2:-1:1
    z_bordeado{zz,1} = padarray(labels_Sc{zz,1},[vec_res2(zz) vec_res2(zz)],nan,'both');
    Z_est =  pattern_0{zz,1}; %padarray(labels_Sc{zz},[1 1],nan,'both');

    TABind0 = find(~isnan(z_bordeado{zz,1}))';

    [TABsub_row, TABsub_col] = ind2sub(size(pattern_0{zz,1}),TABind0);

    clear TABsub_* TABind0

    [TABsub_row{2},TABsub_col{2}]   = find(pattern_0{zz,1}==0);
    [TABsub_row{1},TABsub_col{1}] = find(pattern_0{zz,1}==1);

    TABind{zz,1} = sub2ind(size(Z_est),TABsub_row{1}',TABsub_col{1}');
    n_ind(1) = length(TABind{zz,1});

    TABind{zz,2} = sub2ind(size(Z_est),TABsub_row{2}',TABsub_col{2}');
    n_ind(2) = length(TABind{zz,2}); n_indjj{zz,1} = n_ind;
end




%% 
for zz = j2:-1:2
    patt = ones(1,2);
    D = eye(size(z_bordeado{zz,1}));
    A = kron((D),patt);
    B = kron(patt,A);
    DD = kron(D,B);
    % ----------------------

    for ii = 1:2
        D2 = DD;
        D2(:,TABind{zz-1,ii}) = 100;
        D3 = DD.*D2;
        D3(D3~=100)=0;
        D3(D3==100)=1;

        DD_f{zz,ii} = sparse(D3);
        DD_t{zz,ii} = sparse(D3');
    end
    clear D3;
    clear D2;
end
clear DD





end

