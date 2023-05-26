function [ fig ] = imagesc_custom( XX,YY,FF)
% CUSTOM IMAGESC this custom function is meant to be used exaclty as pcolor
% but plotting array values in pixels (like imagesc) and not in vertices (like pcolor).

fig = imagesc( XX(1,:), YY(:,1), FF);
set(gca,'YDir','normal')        
xticks(XX(1,:));
xticklabels(num2str(XX(1,:)));
yticks(YY(:,1));
yticklabels(num2str(YY(:,1)));
end

