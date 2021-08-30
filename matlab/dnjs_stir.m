function [res] = dnjs_stir(n,j,s)
% Compute the dnjs from Ln*Lj = sum_s dnjs Ls with stirling formula
    stirling = @(n,k) sqrt(2*pi).*n.^(n+.5).*exp(-n) .* (n~=0) + (n==0);
    bin_stir = @(n,k) stirling(n)./stirling(n-k)./stirling(k);
    
    res = 0;
    for n1 = 0:n
        for j1 = 0:j
            for s1 = 0:s
                res = res ...
                    +(-1).^(n1+j1+s1)...
                    .* bin_stir(n,n1) .* bin_stir(j,j1) .* bin_stir(s,s1)...
                    ./ stirling(n1)  ./ stirling(j1)  ./ stirling(s1)...
                    .* stirling(n1+j1+s1);
            end
        end
    end
end