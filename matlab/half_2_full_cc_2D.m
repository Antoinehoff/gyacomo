function [ FF_f ] = half_2_full_cc_2D( FF_h )
%half_2_full_cc_2D Retrieve full domain frequency from half domain
%   Meant for complex conjugate symmetric fields

[Nkr,Nkz] = size(FF_h);

if Nkr > Nkz
    FF_f = zeros(max(size(FF_h)));
    FF_f(1:Nkr,1:Nkz) = FF_h;

    for ikr = 1:Nkr
        for ikz = Nkr/2+2:Nkr
            FF_f(ikr,ikz) = FF_f(Nkr-ikr+1,Nkr-ikz+1);
        end
    end
else
    FF_f = zeros(max(size(FF_h)));
    FF_f(1:Nkr,1:Nkz) = FF_h;

    for ikz = 1:Nkz
        for ikr = Nkz/2+2:Nkz
            FF_f(ikr,ikz) = FF_f(Nkz-ikr+1,Nkz-ikz+1);
        end
    end
end
end

