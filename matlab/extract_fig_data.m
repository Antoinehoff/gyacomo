% tw = [500 1000];
% tw = [1000 1500];
% tw = [2000 3000];
% tw = [3000 4000];
% tw = [4000 4500];
% tw = [4500 5000];
tw = [00 40000];

fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
X_ = [];
Y_ = [];
for i = 1:numel(dataObjs)
    X_ = [X_ dataObjs(i).XData];
    Y_ = [Y_ dataObjs(i).YData];
end
% n0 = 1;
% n1 = numel(X_);
[~,n0] = min(abs(X_-tw(1)));
[~,n1] = min(abs(X_-tw(2)));
mvm = @(x) movmean(x,1);
shift = 0;%X_(n0);
% shift = 0;
skip = 1;

figure;
plot(mvm(X_(n0:skip:n1)-shift),mvm(Y_(n0:skip:n1))); hold on;

% t0 = ceil(numel(X_)*0.2); t1 = numel(X_);
avg= mean(Y_(n0:n1)); dev = std(Y_(n0:n1));
  disp(['AVG =',sprintf('%4.4f',avg),'+-',sprintf('%4.4f',dev)]);
% 
% n1 = n0+1; n2 = min(n1 + 50000,numel(Y_));
% avg_ = mean(Y_(n1:n2));
% std_ = std(Y_(n1:n2));
% title(['avg = ',num2str(avg_),' std = ', num2str(std_)])
% plot([X_(n1),X_(n2)],[1 1]*avg_,'--k')

% ylim([0,avg_*3]);
% sum(Y_(2:end)./X_(2:end).^(3/2))
%%

if 0
    %%
   m1 = mean(y1); m2 = mean(y2);
   s1 = std(y1);  s2 = std(y2);
   delta = abs(m2-m1)/min(s1,s2);
   delta
end


%     xlim([-3,3]); ylim([0,4.5]);