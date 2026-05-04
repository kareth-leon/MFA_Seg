function [ze,ae] = regr_cp(VV,j1,j2,nj,wtype,urs);
if nargin<6; urs=[1 size(y,1)]; end
if numel(urs)==1; urs=[1 urs]; end

%keyboard
jj=j1:j2;
J=length(jj);
njj=nj(jj);
yjj=VV(:,jj);
wvarjj = ones(1,J);

if     wtype==1
    wvarjj = 1./njj;
    %elseif wtype ==2
    %   wvarjj = 1./njj; wvarjj(end)=Inf;
end

% use weighted regression formula in all cases
S0 = sum(1./wvarjj) ;
S1 = sum(jj./wvarjj) ;
S2 = sum(jj.^2./wvarjj) ;
wjj = (S0 * jj - S1) ./ wvarjj / (S0*S2-S1*S1) ;
vjj = (S2 - S1 * jj) ./ wvarjj / (S0*S2-S1*S1) ;


%  Estimate  zeta
%keyboard
%ze  = sum(wjj * yjj' )/log(2) ;       % zeta is just the slope, unbiased regardless of the weights
ze  = (wjj * yjj')/log(2) ;
ae     = (vjj * yjj' ) ;       % intercept  'a'

try ze=reshape(ze,urs); end
try ae=reshape(ae,urs); end