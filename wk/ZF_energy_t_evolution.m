sz_ = size(data.PHI);
iz  = sz_(3)/2+1;
Ezf = sum(abs(data.PHI(:,:,:,:)),2);
Ezf = squeeze(Ezf(1,:,13,:));
Enz = sum(abs(data.PHI(:,:,:,:)),2);
Enz = squeeze(sum(Enz(2:end,:,13,:),1));

figure; 
plot(data.Ts3D,Ezf,'DisplayName','Zonal energy'); hold on;
plot(data.Ts3D,Enz,'DisplayName','Non-zonal energy');
