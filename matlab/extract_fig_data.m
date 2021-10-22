fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;

X_ = dataObjs(1).XData;
Y_ = dataObjs(1).YData;

figure;
plot(X_,Y_);
plot(X_(9000:12000)-X_(8000),Y_(9000:12000));