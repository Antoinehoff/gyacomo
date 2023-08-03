function [w,e,t] = compute_growth_rates(field,time)
% compute_growth_rates compute the linear growth rates using the amplitude
% ratio method
% see B.J. Frei et al. flux-tube paper
% Input: (k1,k2,z,t) field
% Output: w(k1,k2,t) the growth and frequencies w = gamma + i * omega
%         e(k1,k2) an error estimate of the convergence of w
sz = size(field);
N1 = sz(1); N2 = sz(2); Nz = sz(3); Nt = sz(4);


w = zeros(N1,N2,Nt-1);
e = zeros(N1,N2);
for i1 = 1:N1
    for i2 = 1:N2
    to_measure = reshape(field(i1,i2,:,:),Nz,Nt);
        % Amplitude ratio method for determining the growth rates and the
        % frequencies
        for it = 2:Nt
            phi_n   = to_measure(:,it); 
            phi_nm1 = to_measure(:,it-1);
            dt      = time(it)-time(it-1);
            ZS      = sum(phi_nm1,1); %sum over 
            wl          = log(phi_n./phi_nm1)/dt;
            w(i1,i2,it) = squeeze(sum(wl.*phi_nm1,1)./ZS);
        end
        % error estimation
        wavg = mean(w(i1,i2,ceil(Nt/2):end));
        wend = w(i1,i2,end);
        e(i1,i2) = abs(wend-wavg)/abs(wavg);
    end
end

t = time(2:Nt);

end