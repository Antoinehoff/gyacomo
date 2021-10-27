%%

cd ..
system('make');
cd wk

shearless_linear_fluxtube


cd ../../molix
system('make');
system('./bin/molix');
cd ../HeLaZ/wk
%%
time_2_plot = 5.0;
[Y_,X_]=molix_plot_phi('../../molix/field.dat.h5',time_2_plot);

plot_phi_ballooning; hold on
plot(X_/pi,real(Y_),'ob');
plot(X_/pi,imag(Y_),'or');
plot(X_/pi,abs(Y_) ,'ok');