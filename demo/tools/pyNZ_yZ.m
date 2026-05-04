function [yNZ,yZ] = pyNZ_yZ(y,ind0)
% Non-zero frequency
yNZ = y; 
yNZ(:,ind0) = [];
% Separate zero frequency and non zero ones
yZ = y(:,ind0);
end

