function [ ] = plot_gbms_ballooning(resfile)
% perform ballooning transformation from phi.dat.h5 file

% read data and attributes
coordkx = h5read(resfile,'/data/var2d/phi/coordkx');
Nkx = (length(coordkx)-1)/2;
coordz = h5read(resfile,'/data/var2d/phi/coordz');
Nz = length(coordz);
coordtime =h5read(resfile,'/data/var2d/time');
Nt = length(coordtime) ;

iframe = Nt -1;
dataset = ['/data/var2d/phi/',num2str(iframe,'%06d')];
phi = h5read(resfile,dataset);
try
dataset = ['/data/var2d/psi/',num2str(iframe,'%06d')];
psi = h5read(resfile,dataset);
catch
    psi = 0;
end
% Apply baollooning tranform
dims = size(phi.real);
phib_real = zeros(  dims(1)*Nz  ,1);
phib_imag= phib_real;
psib_real= phib_real;
psib_imag= phib_real;
b_angle = phib_real;

midpoint = floor((dims(1)*Nz )/2)+1;

for ip =1: dims(1) 
    start_ =  (ip -1)*Nz +1;
    end_ =  ip*Nz;
    phib_real(start_:end_)  = phi.real(ip,:);  
    phib_imag(start_:end_)  = phi.imaginary(ip,:); 
    try
    psib_real(start_:end_)  = psi.real(ip,:);  
    psib_imag(start_:end_)  = psi.imaginary(ip,:); 
    catch
    psib_real(start_:end_)  = 0;  
    psib_imag(start_:end_)  = 0;
    end
end
    
% Define ballooning angle
idx = -Nkx:1:Nkx;
for ip =1: dims(1)
    for iz=1:Nz
        ii = Nz*(ip -1) + iz;
       b_angle(ii) =coordz(iz) + 2*pi*idx(ip);
    end
end 

% normalize real and imaginary parts at chi =0
[~,idxLFS] = min(abs(b_angle -0));
phib = phib_real + 1i*phib_imag;
psib = psib_real + 1i*psib_imag;
% normalize to the outer mid-plane
norm = (phib(idxLFS));
phib = phib(:)/norm;
psib = psib(:)/norm;
figure;
subplot(1,2,1);
plot(b_angle/pi,real(phib),'b'); hold on;
plot(b_angle/pi,imag(phib),'r'); hold on;
plot(b_angle/pi, abs(phib),'k'); hold on;
xlabel('$\chi/\pi$'); ylabel('$\phi(\chi)/|\phi(0)|$');
 title('GBMS');
subplot(1,2,2)
plot(b_angle/pi,real(psib),'b'); hold on;
plot(b_angle/pi,imag(psib),'r'); hold on;
plot(b_angle/pi, abs(psib),'k'); hold on;
xlabel('$\chi/\pi$'); ylabel('$\psi(\chi)/|\phi(0)|$');
end
