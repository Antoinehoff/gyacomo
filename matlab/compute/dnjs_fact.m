function [res] = dnjs_fact(n,j,s)
% Compute the dnjs from Ln*Lj = sum_s dnjs Ls with factorials only
    binomial = @(n,k) factorial(n)./factorial(n-k)./factorial(k);
    res = 0;
    for n1 = 0:n
        for j1 = 0:j
            for s1 = 0:s
                res = res ...
                    +(-1).^(n1+j1+s1)...
                    .* binomial(n,n1) .* binomial(j,j1) .* binomial(s,s1)...
                    ./ factorial(n1)  ./ factorial(j1)  ./ factorial(s1)...
                    .* factorial(n1+j1+s1);
            end
        end
    end
end