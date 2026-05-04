function NormA = NormA(x,A)
% Given x a nxm array and A a mxm matrix, returns a nx1 column vector
% where the i-th row contains x(i,:)*A*x(i,:)'
% m = size(A,1);
% if size(x,1)==m
%     x = transpose(x);
% end
% tic
% NormA = (diag(x*A*x'));
% toc
%   keyboard
%keyboard
if nargin >1
    if numel(A)~= length(A)
        m = size(A,1);
        if size(x,2)==m
            x = transpose(x);
        end
        NormA = real(sum(transpose(x').*(A*x)))';
    else
        NormA = (sum(bsxfun(@times,abs(x).^2,A'),2));
    end    
else
    NormA = sum(abs(x).^2,2);  
end
end