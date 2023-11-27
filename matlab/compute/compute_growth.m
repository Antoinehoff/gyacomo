function [gamma,err] = compute_growth(t,y)
gamma = zeros(size(y));
    for it = 2:numel(t)
        y_n   = y(it); 
        y_nm1 = y(it-1); 
        dt    = t(it)-t(it-1);
        ZS    = sum(y_nm1,1);
        wl    = log(y_n./y_nm1)/dt;
        gamma(it) = squeeze(sum(wl.*y_nm1,1)./ZS);
    end
    % error estimation
    n          = 5;
    [gamma, ~, err] = sliceAverage(gamma, n);
    err = mean(err);
end