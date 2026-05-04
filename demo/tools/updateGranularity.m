function betaNew = updateGranularity(z,z_k,pmlbl,loglike_border,vec_res2,Njj,iter2)
%========================================
% Author: Kareth LEON
% Function: update the granularity function
% 
%========================================


betaNew =zeros(size(pmlbl.beta_current));

newZ = cell(pmlbl.j2,1);
ones_like = cell(pmlbl.j2,1);

for ii = 1:pmlbl.j2
    newZ{ii,1} = mod(z{ii}+1,pmlbl.K) + 1 ;
    ones_like{ii,1} = ones(size(loglike_border{ii}));
end

%% random generation of labels
for tt = 1:pmlbl.it
    newZ = updateChessboard_MATLAB_SYS(newZ,ones_like,pmlbl);
end
idxZ = z_k; 
idxNewZ = createInd(newZ,pmlbl,vec_res2);


%% 
[n,m] = size(pmlbl.beta_current);
A = pottsLogLikelihood_MR(idxZ,ones(n,m),pmlbl);
B = pottsLogLikelihood_MR(idxNewZ,ones(n,m),pmlbl);

betaOld = pmlbl.beta_current;

%% CONSTANT
eta_j = 10*((sum(prod(Njj'))).^-1)*(iter2-1)^(-3/4);

%% ESTIMATION
for zz = 1:pmlbl.j2
    % SPATIAL BETA
    betaNew(1,zz) = betaOld(1,zz) + eta_j*(A(3,1)-B(3,1));
    betaNew(3,zz) = betaOld(3,zz) + eta_j*(A(3,1)-B(3,1));
end



if (sum(find(betaNew > 10)))
    betaNew(find(betaNew> 10))=10;
elseif(sum(betaNew(:)<0))
    betaNew(find(betaNew<0))=0;
end
end



%%
function aux_p = pottsLogLikelihood_MR(z,beta,pmlbl)

% compute derivatives
for zz=1:pmlbl.j2
    aux = z{zz};
    vd{zz}=diff(aux);         % compute vertical differences
    hd{zz}=diff(aux,1,2);     % compute horizontal differences

    if zz < pmlbl.j2 % compute scale differences
        for kk = 1:pmlbl.K
            z_next_replicate = repelem(z{zz+1}(:,:,kk),2,2);
            aux_d = cat(3,aux(:,:,kk),z_next_replicate);
            dd{zz}(:,:,kk) =  diff(aux_d,1,3);
        end
    else
        dd{zz} = 0.0001;
    end
end

% sum up changes
auxx = 0;
for zz=1:pmlbl.j2
    aux_p(1,zz) = beta(1,zz)*sum(~vd{zz}(:)) + beta(1,zz)*sum(~hd{zz}(:));
    auxx = auxx + aux_p(1,zz) +  beta(3,zz)*(sum(~dd{zz}(:)));
end
aux_p(3,1) = auxx;
end



function Z_est_cut = createInd(oldZ,pmlbl,vec_res2)

% function to create K version of the labels
for zz = 1:pmlbl.j2
    for kk = 1:pmlbl.K
        Z_est_cut{zz,1}(:,:,kk) = oldZ{zz,1}(vec_res2(zz)+1:end-vec_res2(zz),vec_res2(zz)+1:end-vec_res2(zz))==kk;
    end
end
end

