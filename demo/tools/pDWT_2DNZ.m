function [Nj,ind,VV,mV,Lx] = pDWT_2DNZ(im3D,j1,j2,Nwt,eta,gamint)
%==========================================================================
% function from toolbox :
% https://www.irit.fr/~Herwig.Wendt/CIMPA2017.html
%==========================================================================

if nargin < 6
    gamint = 0;
end
M = size(im3D,3);

for t=1:M

    tmp=im3D(:,:,t);
    if sum(~isfinite(tmp(:)))
        for j=j1:j2
            leaders(j).matrix=inf(size(im3D(1:2^j:end,1:2^j:end,t)));
        end
    else
        % Wavelet transform
        [~, leaders, ~] = DxLx2d_bayes_dwt(im3D(:,:,t),Nwt,j2,gamint);%Nj = sqrt(nj.L);
        if length(leaders)<j2 && ~sum(isfinite(leaders(end).allmatrix(:)))
            for j=length(leaders)+1:j2; leaders(j).allmatrix=leaders(j-1).allmatrix(1:2:end,1:2:end); end
        end
    end

    if length(leaders)<j2
        fprintf('\n\n ERROR:\n There are less than j2 scales available.\n Reduce j2 and restart.\n\n');
        error(' ');
    end

    for j=j1:j2
        lx{j}   = log(leaders(j).matrix);
        Lx{j}{t} = lx{j};
        VV(t,j) = var(lx{j}(:));
        mV(t,j) = mean(lx{j}(:));
        Nj(j,:) = size(lx{j});
    end

    if t==1;[~,ind] = fft_freqNZ(j1,j2,Nj,eta);end


end
end










