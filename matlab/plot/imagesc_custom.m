function [ fig ] = imagesc_custom( XX,YY,FF)
% CUSTOM IMAGESC this custom function is meant to be used exaclty as pcolor
% but plotting array values in pixels (like imagesc) and not in vertices (like pcolor).
x = reshape(XX(1,:),[numel(XX(1,:)), 1]);
y = reshape(YY(:,1),[numel(YY(:,1)), 1]);
fig = imagesc( x, y, FF);
set(gca,'YDir','normal')        
xticks(x);
xticklabels(num2str(x));
yticks(y);
yticklabels(num2str(y));
end

