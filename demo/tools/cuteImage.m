function [im0,C,u0,b0] = cuteImage(im,P,dP)
%keyboard

[M,N,~] = size(im);
rr = P(1); cc = P(2);
xx = dP(1); yy = dP(2);

[XX,YY]=ndgrid( 1:(rr-xx):(M-(rr-1)) , 1:(cc-yy):(N-(cc-1)));

IDX1 = bsxfun(@plus, (1:rr)', ((1:cc)-1)*M);

% offset in indices for each block
offset = (XX(:)-1) + (YY(:)-1)*M;
% indices of each block, the block is indexed with the 3rd dimension
IDX = bsxfun(@plus, IDX1, reshape(offset, [1 1 numel(offset)]));
% convert to cell of blocks

im0 =  im(IDX);
u0 = size(XX);
if prod(u0)>1
    C = reshape(mat2cell(im(IDX), rr, cc, ones(1, numel(XX))),u0);
else
    %keyboard
    warning off
    C = mat2cell(im(IDX), rr, cc, ones(1, numel(XX)));
    warning on
end

b0 = [max(XX(:)) max(YY(:))]+P-1;
end
