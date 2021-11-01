fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;

X_ = dataObjs(1).XData;
Y_ = dataObjs(1).YData;

figure;
plot(X_,Y_);
plot(X_-X_(1),Y_);