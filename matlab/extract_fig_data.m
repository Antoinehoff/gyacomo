fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;

X_ = dataObjs(1).XData;
Y_ = dataObjs(1).YData;
n0 = 0001;
figure;
plot(X_(n0:end),Y_(n0:end));
plot(X_(n0:end)-X_(n0),Y_(n0:end));