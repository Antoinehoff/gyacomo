function [ result ] = dnjs( n, j, s )
% Compute the dnjs from Ln*Lj = sum_s dnjs Ls with Laguerre coeffs
% sort in order to compute only once the laguerre coeff
    Coeffs = sort([n,j,s]); 
% last element of Coeffs is the larger one
    L3 = flip(LaguerrePoly(Coeffs(end)));
    L2 = flip(LaguerrePoly(Coeffs(end-1)));
    L1 = flip(LaguerrePoly(Coeffs(end-2)));

% build a factorial array to compute once every needed factorials
    Factar    = zeros(n+j+s+1,1);
    Factar(1) = 1.0;
    f_ = 1;
    for ii_ = 2:(n+j+s)+1
        f_ = (ii_-1) * f_;
        Factar(ii_) = f_;
    end
    
    result = 0;
    for il3 = 1:numel(L3)
        for il2 = 1:numel(L2)
            for il1 = 1:numel(L1)
                result = result + Factar(il1 + il2 + il3 -2)...
                        *L1(il1) * L2(il2) * L3(il3); 
            end
        end
    end
end