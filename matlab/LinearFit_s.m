function [gamma,fit] = LinearFit_s(time,Na00abs)
  % LinearFit computes the growth rate and frequency from the time evolution of Napj
  % - adapted from MOLI (B.J. Frei)
  %

    % ... amplitude ratio method

    % We compute the mean of the growth rate over a time window
%     Trun = time(end);
%     if tstart == -1
%         [~,Tstart_ind] = min(abs(time - 0.5*Trun));
%     else
%         [~,Tstart_ind] = min(abs(time - max(tstart,time(1))));
%     end
%     
%     if tend == -1
%         Tend_ind       = numel(time)-1;
%     else
%         [~,Tend_ind]   = min(abs(time - min(tend,Trun)));
%     end
%         


    Na00absshifted = circshift(Na00abs,-1);    % ... shift by -1 the time position 
    gammaoft = log(Na00absshifted(1:end-1)./Na00abs(1:end-1))./squeeze(diff(time')); % ... evaluate growth rate

    % Get gamma
%     gamma = mean(gammaoft(Tstart_ind:Tend_ind-1)); % ... take the mean of gamma over the time window
    gamma = mean(gammaoft); % ... take the mean of gamma over the time window

    % Return gamma(t) for amplitude ratio method
    fit.gammaoft = gammaoft;

end % ... end function
