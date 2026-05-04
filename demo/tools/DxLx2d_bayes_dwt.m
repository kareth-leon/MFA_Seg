function [coef, leaders, nj] = DxLx2d_bayes_dwt(data, Nwt,JM, gamint,symm)
%==========================================================================
%
%  [coef, leaders, nj] = DxLx2d(data, Nwt [,gamint,sym])
%
% compute Dyadic wavelet and leader coefficients  
%
% gamint -- fractional integration of order gamint is applied
% sym=0  -- use Daubechies Wavelet of order Nwt (default).
% sym=1  -- use Daubechies Symmetrized Wavelet of order Nwt.
%
%
%   USAGE EXAMPLES :
%
% N=2^9;
% Nwt=5;
% signal=randn(N,N);
% [coef, leaders, nj] =  DxLx2d(signal, Nwt);
%
%%%
% SR, ens-Lyon, 11/2013
%
%%%%%%%%%%
% MODIF 
%leaders(j).vector_clean : vecteur des leaders avec retrait des leaders
%infinis 
% SC 03/06/14
%%%%%%%%%%
%==========================================================================

Norm=1;

%-- paramater for  fractional integration
if nargin<4; gamint=0; end;
if isempty(gamint); gamint=0; end;

if nargin<5; symm=0; end;
%-- Initialize the wavelet filters
n = length(data) ;         % data length

if symm == 0 % Daubechies Wavelet
    h = rlistcoefdaub(Nwt) ;   % filter
    nl = length(h) ;           % length of filter, store to manage edge effect later
    gg1 = -1*(-1).^(1:nl).*h ; % wavelet filter
    hh1=fliplr(h);             % scaling filter
    % parameter for the centering of the wavelet
    x0=2; 
    x0Appro=2*Nwt; 

else % Daubechies Symmetrized Wavelet
    Nwt=abs(Nwt);
    h = rlistcoefdaub(Nwt) ;   % filter
    nl = length(h) ;           % length of filter, store to manage edge effect later
    gg1c = -1*(-1).^(1:nl).*h ; % wavelet filter
    hh1c = fliplr(h);             % scaling filter
    tmp = conv(h,hh1c)/sqrt(2);
    nu = -2*Nwt+1;
    [v,nv] = tildeurflippeur(tmp,nu);
    
    [hh1,nh1] = flippeur(tmp,nu) ; nl=length(hh1);
    [gg1,ng1] = flippeur(v,nv) ;   

    
    % parameter for the centering of the wavelet
    x0=2*Nwt; 
    x0Appro=2*Nwt; 

end
%--- Predict the max # of octaves available given Nwt, and take the min with
nbvoies= fix( log2(length(data)) );
nbvoies = min( fix(log2(n/(nl+3)))  , nbvoies); %   safer, casadestime having problems

if nargin>2;
    nbvoies =min(nbvoies, JM);
end

coef=struct([]);
%--- Compute the WT, calculate statistics
LL=data;
sidata=size(data);

for j=1:nbvoies         % Loop Scales
    
    njtemp = size(LL) ;
    %-- border effect
    fp=nl; % index of first good value
    lp=njtemp; % index of last good value
    
    %-- OH convolution and subsampling
    OH=conv2(LL,gg1); OH(isnan(OH))=Inf;
    OH(:,1:fp-1)=Inf;
    OH(:,lp(2)+1:end)=Inf;
    OH=OH(:,(1:2:njtemp(2))+x0-1);
    
    %-- HH convolution and subsampling
    HH=conv2(OH,gg1');HH(isnan(HH))=Inf;
    HH(1:fp-1,:)=Inf;
    HH(lp(1)+1:end,:)=Inf;
    HH=HH((1:2:njtemp(1))+x0-1,:);

    %-- LH convolution and subsampling
    LH=conv2(OH,hh1');LH(isnan(LH))=Inf;
    LH(1:fp-1,:)=Inf;
    LH(lp(1)+1:end,:)=Inf;
    LH=LH((1:2:njtemp(1))+x0Appro-1,:);
    clear OH
    
    %-- OL convolution and subsampling
    OL=conv2(LL,hh1);OL(isnan(OL))=Inf;
    OL(:,1:fp-1)=Inf;
    OL(:,lp(2)+1:end)=Inf;
    OL=OL(:,(1:2:njtemp(2))+x0Appro-1);

    %-- HL convolution and subsampling
    HL=conv2(OL,gg1');HL(isnan(HL))=Inf;
    HL(1:fp-1,:)=Inf;
    HL(lp(1)+1:end,:)=Inf;
    HL=HL((1:2:njtemp(1))+x0-1,:);

    %-- LL convolution and subsampling
    LL = conv2(OL,hh1');LL(isnan(LL))=Inf;
    LL(1:fp-1,:) = Inf;
    LL(lp(1)+1:end,:) = Inf;

    LL = LL((1:2:njtemp(1))+x0Appro-1,:);
    clear OL
    
    %-- passage Norme L1
    ALH=abs(LH)/2^(j/Norm);
    AHL=abs(HL)/2^(j/Norm);
    AHH=abs(HH)/2^(j/Norm);
    ALL=abs(LL)/2^(j/Norm); % SC 2013
    
        
    %-- fractional integration by gamma
    coef(j).gamma=gamint;
    ALH=ALH*2^(gamint*j);
    AHL=AHL*2^(gamint*j);
    AHH=AHH*2^(gamint*j);
    
    coef(j).allx=ALH;
    coef(j).ally=AHL;
    coef(j).allxy=AHH;

    %-- get position of coefs
    lesx=1:2^j:sidata(2);
    lesy=1:2^j:sidata(1);
    
    [i1,j1]=find(isfinite(ALH));
    [i2,j2]=find(isfinite(AHL));
    [i3,j3]=find(isfinite(AHH));   
    ii1=max([min(i1) min(i2) min(i3)]);
    ii2=min([max(i1) max(i2) max(i3)]);
    jj1=max([min(j1) min(j2) min(j3)]);
    jj2=min([max(j1) max(j2) max(j3)]);
    
    coef(j).xpos=lesx(jj1:jj2);
    coef(j).ypos=lesy(ii1:ii2);
    coef(j).x=ALH(ii1:ii2,jj1:jj2);
    coef(j).y=AHL(ii1:ii2,jj1:jj2);
    coef(j).xy=AHH(ii1:ii2,jj1:jj2);
    coef(j).approx=ALL(ii1:ii2,jj1:jj2);
     
    coef(j).supcoef=max([coef(j).x(:)' coef(j).y(:)' coef(j).xy(:)']);
    
    nj.W(j)=3*prod(size(coef(j).x));

    if j ==1
        %-- compute and store leaders sans voisin
        leaders(j).sans_voisin.value.allx = coef(j).allx;
        leaders(j).sans_voisin.value.ally = coef(j).ally;
        leaders(j).sans_voisin.value.allxy = coef(j).allxy;        
    else
        nc=floor(size(leaders(j-1).sans_voisin.value.allx)/2);
        %-- get max at smaller scales
        % x
        clear temp*
        temp1x=coef(j).allx(1:nc(1),1:nc(2));
        temp2x=leaders(j-1).sans_voisin.value.allx(1:2:2*nc(1),1:2:2*nc(2)); 
        temp3x=leaders(j-1).sans_voisin.value.allx(2:2:2*nc(1),1:2:2*nc(2));
        temp4x=leaders(j-1).sans_voisin.value.allx(1:2:2*nc(1),2:2:2*nc(2));
        temp5x=leaders(j-1).sans_voisin.value.allx(2:2:2*nc(1),2:2:2*nc(2));
        tempx(1,:,:)=temp1x; tempx(2,:,:)=temp2x; tempx(3,:,:)=temp3x; tempx(4,:,:)=temp4x; tempx(5,:,:)=temp5x;
        leaders(j).sans_voisin.value.allx=squeeze(max(tempx));

        % y
        clear temp*
        temp1x=coef(j).ally(1:nc(1),1:nc(2));
        temp2x=leaders(j-1).sans_voisin.value.ally(1:2:2*nc(1),1:2:2*nc(2)); 
        temp3x=leaders(j-1).sans_voisin.value.ally(2:2:2*nc(1),1:2:2*nc(2));
        temp4x=leaders(j-1).sans_voisin.value.ally(1:2:2*nc(1),2:2:2*nc(2));
        temp5x=leaders(j-1).sans_voisin.value.ally(2:2:2*nc(1),2:2:2*nc(2));
        tempx(1,:,:)=temp1x; tempx(2,:,:)=temp2x; tempx(3,:,:)=temp3x; tempx(4,:,:)=temp4x; tempx(5,:,:)=temp5x;
        leaders(j).sans_voisin.value.ally=squeeze(max(tempx));
        
        % xy
        clear temp*
        temp1x=coef(j).allxy(1:nc(1),1:nc(2));
        temp2x=leaders(j-1).sans_voisin.value.allxy(1:2:2*nc(1),1:2:2*nc(2)); 
        temp3x=leaders(j-1).sans_voisin.value.allxy(2:2:2*nc(1),1:2:2*nc(2));
        temp4x=leaders(j-1).sans_voisin.value.allxy(1:2:2*nc(1),2:2:2*nc(2));
        temp5x=leaders(j-1).sans_voisin.value.allxy(2:2:2*nc(1),2:2:2*nc(2));
        tempx(1,:,:)=temp1x; tempx(2,:,:)=temp2x; tempx(3,:,:)=temp3x; tempx(4,:,:)=temp4x; tempx(5,:,:)=temp5x;
        leaders(j).sans_voisin.value.allxy=squeeze(max(tempx));
        
        
        
    end
   
    %-- on prend le max sur les 8 voisins i.e. 9 coeffs
    % x
    clear lsx temp*;
    six=size(leaders(j).sans_voisin.value.allx);
    lsx = zeros(2+six(1),2+six(2));
    tempx = zeros(six(1),six(2),9);
    lsx(2:end-1,2:end-1) = leaders(j).sans_voisin.value.allx;
    tempx(:,:,1) = lsx(1:end-2,1:end-2);
    tempx(:,:,2) = lsx(1:end-2,2:end-1);
    tempx(:,:,3) = lsx(1:end-2,3:end);
    tempx(:,:,4) = lsx(2:end-1,1:end-2);
    tempx(:,:,5) = lsx(2:end-1,2:end-1);
    tempx(:,:,6) = lsx(2:end-1,3:end);
    tempx(:,:,7) = lsx(3:end,1:end-2);
    tempx(:,:,8) = lsx(3:end,2:end-1);
    tempx(:,:,9) = lsx(3:end,3:end);
    leaders(j).value.allx = max(tempx,[],3);
    
    % y
    clear lsx temp*;
    six=size(leaders(j).sans_voisin.value.ally);
    lsx = zeros(2+six(1),2+six(2));
    tempx = zeros(six(1),six(2),9);
    lsx(2:end-1,2:end-1) = leaders(j).sans_voisin.value.ally;
    tempx(:,:,1) = lsx(1:end-2,1:end-2);
    tempx(:,:,2) = lsx(1:end-2,2:end-1);
    tempx(:,:,3) = lsx(1:end-2,3:end);
    tempx(:,:,4) = lsx(2:end-1,1:end-2);
    tempx(:,:,5) = lsx(2:end-1,2:end-1);
    tempx(:,:,6) = lsx(2:end-1,3:end);
    tempx(:,:,7) = lsx(3:end,1:end-2);
    tempx(:,:,8) = lsx(3:end,2:end-1);
    tempx(:,:,9) = lsx(3:end,3:end);
    leaders(j).value.ally = max(tempx,[],3);
    
    
    % xy
    clear lsx temp*;
    six=size(leaders(j).sans_voisin.value.allxy);
    lsx = zeros(2+six(1),2+six(2));
    tempx = zeros(six(1),six(2),9);
    lsx(2:end-1,2:end-1) = leaders(j).sans_voisin.value.allxy;
    tempx(:,:,1) = lsx(1:end-2,1:end-2);
    tempx(:,:,2) = lsx(1:end-2,2:end-1);
    tempx(:,:,3) = lsx(1:end-2,3:end);
    tempx(:,:,4) = lsx(2:end-1,1:end-2);
    tempx(:,:,5) = lsx(2:end-1,2:end-1);
    tempx(:,:,6) = lsx(2:end-1,3:end);
    tempx(:,:,7) = lsx(3:end,1:end-2);
    tempx(:,:,8) = lsx(3:end,2:end-1);
    tempx(:,:,9) = lsx(3:end,3:end);
    leaders(j).value.allxy = max(tempx,[],3);
    
   
    % get the position of leaders
    [i1,j1]=find(isfinite(leaders(j).value.allx));
    [i2,j2]=find(isfinite(leaders(j).value.ally));
    [i3,j3]=find(isfinite(leaders(j).value.allxy));
    
    
    ii1=max([min(i1) min(i2) min(i3)]);
    ii2=min([max(i1) max(i2) max(i3)]);
    jj1=max([min(j1) min(j2) min(j3)]);
    jj2=min([max(j1) max(j2) max(j3)]);

        
    leaders(j).value.xpos=lesx(jj1:jj2);
    leaders(j).value.ypos=lesy(ii1:ii2);
    leaders(j).value.xid=(jj1:jj2);
    leaders(j).value.yid=(ii1:ii2);
    leaders(j).value.x=leaders(j).value.allx(ii1:ii2,jj1:jj2);
    leaders(j).value.y=leaders(j).value.ally(ii1:ii2,jj1:jj2);
    leaders(j).value.xy=leaders(j).value.allxy(ii1:ii2,jj1:jj2);        
    
    clear tempmat;
    tempmat(:,:,1)=leaders(j).value.x;
    tempmat(:,:,2)=leaders(j).value.y;
    tempmat(:,:,3)=leaders(j).value.xy;
    leaders(j).matrix=max(tempmat,[],3);

    clear tempmat;
    tempmat(:,:,1)=leaders(j).value.allx;
    tempmat(:,:,2)=leaders(j).value.ally;
    tempmat(:,:,3)=leaders(j).value.allxy;
    leaders(j).allmatrix=max(tempmat,[],3);
    % SC 03/06/14
    
    leaders(j).gamma=gamint;


    nj.L(j) = numel(leaders(j).value.x);


    
end

%in case of last scale empty
if isempty(leaders(end).value.x)
    coef=coef(1:end-1);
    leaders=leaders(1:end-1);
    nj.W=nj.W(1:end-1);
    nj.L=nj.L(1:end-1);
end
