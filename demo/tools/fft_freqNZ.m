function [w,ind] = fft_freqNZ(j1,j2,Nj,eta)
% Compute frequencies used in fft
% works in 1D and 2D
% keyboard

for j=j1:j2

    fs = 2*pi;
    df1 = fs/Nj(j,1); 
    wk1 = (0:df1:(fs-df1)) - (fs-mod(Nj(j,1),2)*df1)/2; 
    n01 = floor(Nj(j,1)/2)+1;
    
    
    df2 = fs/Nj(j,2); 
    wk2 = (0:df2:(fs-df2)) - (fs-mod(Nj(j,2),2)*df2)/2; 
    n02 = floor(Nj(j,2)/2)+1;
    
    [wx,wy] = meshgrid(wk1,wk2);wr = sqrt(wx.^2+wy.^2);

    ind{j} = ((wr<=sqrt(eta)*pi)&(wr>0.001)&(wx>=0)&(wy>=0))|((wr<=sqrt(eta)*pi)&(wr>0.001)&(wx>0.0001)&(wy<0));
    ind{j}(n02,n01) = 1;


    %keyboard
    %    ind{j} = ((wr<sqrt(eta)*pi)&(wy>0));  ind{j}(n02,n01)=1;
    % ind{j} = ((wr<sqrt(eta)*pi)&(wx>=0)&(wy>0.1));

    w.x{j} = wx(ind{j}); % First composant
    w.y{j} = wy(ind{j}); % Second composant
    w.r{j} = wr(ind{j}); % Radius in 2D / absolute value in 1D
end


end

