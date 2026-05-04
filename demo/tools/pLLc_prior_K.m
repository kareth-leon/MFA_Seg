function [mm,Fi,Gi,iFi,iGi,pL] = pLLc_prior_K(LLc,j1,j2,Nj,ind0, U, K)
%==========================================================================

%==========================================================================

Fi = cell(K,1);
iFi = cell(K,1);
Gi = cell(K,1);
iGi = cell(K,1);
%.-----------------------------------------------------------------------
for kk = 1:K
    auxLLc = LLc{kk};

    indM = find(-auxLLc(1,:)<0);
    m  = min(-auxLLc(2,indM)./(-auxLLc(1,indM)))-10^(-3.5);
    mm{kk}  = m;

    % Non-zero frequency
    LLcNZ = auxLLc;
    LLcNZ(:,ind0{kk}) = [];

    fi  = -LLcNZ(1,:);
    gi = LLcNZ(2,:);

    Fi{kk}  = (m * fi + gi )/m;
    Gi{kk}  = gi / m;

    % Separate zero frequency and non zero ones
    LLcZ = auxLLc(:,ind0{kk});

    %f0 = -LLcZ(1,:);
    %g0 = LLcZ(2,:);

    %F0i  = (m * f0 + g0) / m;
    %G0i  = g0 / m;

    iFi{kk} = 1./Fi{kk};    iGi{kk} = 1./Gi{kk};
end

%==========================================================================

% Prior parameters
pL.param.alpha = 10^-3; pL.param.beta = 10^-3;
pL.param.w0 = zeros(2,1);
pL.param.D = 10^(5)*[1;1];

try pL.EST.Nbi=Nbi; catch pL.EST.Nbi=300; end
try pL.EST.Nmc=Nmc;  catch pL.EST.Nmc=1000; end
try pL.EST.GMRF.a = hpGMRF(1)*ones(1,4); catch pL.EST.GMRF.a = 1.5*ones(1,4); end % Spatial
try pL.EST.GMRF.b = hpGMRF(2)*ones(1,2); catch pL.EST.GMRF.b = 40*ones(1,2); end % Temporal
try pL.EST.GMRF.c = hpGMRF(3)*ones(1,8); catch pL.EST.GMRF.c = 2*ones(1,8); end % 3D anisotropique
end

