function [ny,yNZ_k,M] = fourier_est(I_K,ind0,K)

for kj = 1:K
    [yNZ,yZ] = pyNZ_yZ(I_K{kj}.Y,ind0{kj});
    yNZ_k{kj,1}=yNZ;
    yZ_k{kj,1}=yZ;
end
ny = size(yNZ_k{1},2);      
M = 1; 
