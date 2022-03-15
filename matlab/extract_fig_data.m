fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
X_ = [];
Y_ = [];
for i = 1:numel(dataObjs)
    X_ = [X_ dataObjs(i).XData];
    Y_ = [Y_ dataObjs(i).YData];
end
n0 = 1;
figure;
mvm = @(x) movmean(x,1);
shift = X_(n0);
% shift = 0;
% plot(X_(n0:end),Y_(n0:end));
plot(mvm(X_(n0:end)-shift),mvm(Y_(n0:end))); hold on;

% 
% n1 = n0+1; n2 = min(n1 + 50000,numel(Y_));
% avg_ = mean(Y_(n1:n2));
% std_ = std(Y_(n1:n2));
% title(['avg = ',num2str(avg_),' std = ', num2str(std_)])
% plot([X_(n1),X_(n2)],[1 1]*avg_,'--k')

% ylim([0,avg_*3]);
sum(Y_(2:end)./X_(2:end).^(3/2))
%%

if 0
    %%
   m1 = mean(y1); m2 = mean(y2);
   s1 = std(y1);  s2 = std(y2);
   delta = abs(m2-m1)/min(s1,s2);
   delta
end


%     xlim([-3,3]); ylim([0,4.5]);