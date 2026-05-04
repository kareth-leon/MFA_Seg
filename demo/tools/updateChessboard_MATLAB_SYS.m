function Z_lbl = updateChessboard_MATLAB_SYS(Z_lbl,loglike_border,pmlbl)
%========================================
% Author: Kareth LEON, Abderrahim HALIMI
% Function:
%
%========================================
TABind      = pmlbl.TABind;
n_indjj     = pmlbl.n_indjj;
j2          = pmlbl.j2;


% ------------------------------------------------
% Logic: First all the "black" voxels
% ------ Then, scan and update the "white" voxels
% ------------------------------------------------

% -------------------------------------------------------------------
% ------------------- Chessboard 3D ---------------------------------
%---------------------------------------------------------------------
% WHITE voxels - black/white
%---------------------------------------------------------------------
for zz = 1: 2 : j2  % odd voxels across dim 3
    indid = 1;       % even voxels across dim 1 y 2
    n_ind = n_indjj{zz,1};
    ind3  = TABind{zz,indid};
    Z_lbl{zz,1} = updateVoxel_SYS(Z_lbl,loglike_border,zz,indid,n_ind,ind3,pmlbl);
end

%---------------------------------------------------------------------
% BLACK voxels - black/black
%---------------------------------------------------------------------
for zz = 1:2:j2    % odd voxels across dim 3
    indid = 2;     % odd voxels across dim 1 y 2
    n_ind = n_indjj{zz,1};
    ind3 = TABind{zz,indid};
    Z_lbl{zz,1} = updateVoxel_SYS(Z_lbl,loglike_border,zz,indid,n_ind,ind3,pmlbl);
end

%---------------------------------------------------------------------
% WHITE voxels - white/white
%---------------------------------------------------------------------
for zz = 2:2:j2    % odd voxels across dim 3
    indid   = 1;     % odd voxels across dim 1 y 2
    n_ind   = n_indjj{zz,1};
    ind3    = TABind{zz,indid};
    Z_lbl{zz,1} = updateVoxel_SYS(Z_lbl,loglike_border,zz,indid,n_ind,ind3,pmlbl);
end

%
%---------------------------------------------------------------------
% BLACK voxels - white/black
%---------------------------------------------------------------------
for zz = 2 : 2 : j2  % odd voxels across dim 3
    indid = 2;       % even voxels across dim 1 y 2
    n_ind = n_indjj{zz,1};
    ind3 = TABind{zz,indid};
    Z_lbl{zz,1} = updateVoxel_SYS(Z_lbl,loglike_border,zz,indid,n_ind,ind3,pmlbl);
end
%



end