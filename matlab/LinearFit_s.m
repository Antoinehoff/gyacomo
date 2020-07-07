function [gamma,fit] = LinearFit_s(time,Na00abs)
  % LinearFit computes the growth rate and frequency from the time evolution of Napj
  % - adapted from MOLI (B.J. Frei)
  %

    % ... amplitude ratio method

    % We compute the mean of the growth rate over a time window [0.8*Trun,]
    Trun = time(end);

    lowerbound_timewindow = 0.8*Trun;
    [~,begin_timewindow_ind] = min(abs(time - lowerbound_timewindow));

    Na00absshifted = circshift(Na00abs,-1);    % ... shift by -1 the time position 
    gammaoft = log(Na00absshifted(1:end-1)./Na00abs(1:end-1))./(diff(time)); % ... evaluate growth rate

    % Get gamma
    gamma = mean(gammaoft(end-begin_timewindow_ind:end)); % ... take the mean of gamma over the time window

    % Return gamma(t) for amplitude ratio method
    fit.gammaoft = gammaoft;


    % Return fit
    fit.t_fit_min = lowerbound_timewindow;
    fit.t_fit_max = Trun;

end % ... end function
