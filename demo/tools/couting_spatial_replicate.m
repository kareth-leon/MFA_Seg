function x_out = couting_spatial_replicate(z)
%==========================================================================

%==========================================================================

z_padd = padarray(z,[1 1],nan,'both');
[m,n] = size(z_padd);

j = 1;

for ii = 1:m-2
    i = 1;
    for jj = 1:n-2
        x = z_padd(ii : (ii)+ 2 ,jj : (jj) + 2);
        padre = x(2,2);
        x = x(:);    
        x = x(~isnan(x));
        [a,b] = hist(x,unique(x));
       
        if isscalar(unique(a))
            valor = padre;
        else
            [~,jr] = (max(a));
            valor = b(jr);
        end
        x_out(j:j+1,i:i+1) = repmat(valor,[2,2]);
        i=i+2;

    end
    j=j+2;
end
