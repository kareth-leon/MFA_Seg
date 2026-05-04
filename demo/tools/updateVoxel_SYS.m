function output = updateVoxel_SYS(Z_lbl,loglike_border,zz,indid,n_ind,ind3,pmlbl)
%========================================
% Author: Kareth LEON, Abderrahim HALIMI
% Function:
%
%========================================

DD_t = pmlbl.DD_t;
DD_f = pmlbl.DD_f;
NjjjNUEVA = pmlbl.NjjjNUEVA;


j2 = pmlbl.j2;
K = pmlbl.K;
beta_current = pmlbl.beta_current;



pk_aux = [];
zz_depth = [];

% Search the centers at the third dimension
if zz > 1       % no estoy en el borde
    zz_depth = [zz_depth,zz-1];
end
if zz < j2      % no estoy en el borde
    zz_depth = [zz_depth,zz+1];
end

no_depth = 0;

% iff there is just one scale
if isempty(zz_depth);no_depth = 1;end %just in case, no depth


for kk= 1:K
    sum_delta = 0;
    loglike = loglike_border{zz,1}(:,:,kk);

    % -------------------------------------------------------------
    % Spatial
    % -------------------------------------------------------------
    V           = [ind3-1; ind3+1; ind3-(NjjjNUEVA{zz}(1));ind3+(NjjjNUEVA{zz}(2))];
    sum_delta   = sum_delta + beta_current(1,zz) * sum(kron_index(kk - Z_lbl{zz,1}(V))) ;
    % -------------------------------------------------------------
    % Scales
    % -------------------------------------------------------------
    if no_depth == 0
        long_z = length(zz_depth);
        for z_int = 1 : long_z

            z_d = zz_depth(z_int); % back or front
            
            if zz < z_d
                z_kk        = Z_lbl{z_d,1}(:) == kk;
                valores2    = DD_t{z_d,indid}*z_kk; % salida tiene el tamaño de z_d
                b_mapped    = valores2(ind3)';

                % beta (3,2) --- from 2 to 3
                b_delta_mapped  =  beta_current(3,zz) *b_mapped;

            else
                z_kk        = Z_lbl{z_d,1}(:) == kk;
                producto    = DD_f{zz,indid}*z_kk;
                val1        = max(producto);
                valores     = ((1/val1)*producto); 
                b_mapped    = valores(ind3)';

                b_delta_mapped  =  beta_current(3,z_d) *b_mapped;


            end
            sum_delta = sum_delta + b_delta_mapped;
            clear delta_depth
        end

    end
    % Total sum
    pk = exp( sum_delta + loglike(ind3));
    pk_aux = [pk_aux;pk];
end

output = Z_lbl{zz,1};
% Sampling
output(ind3) = gen_discrete_parallel(1:K,pk_aux,1,n_ind(indid));

end