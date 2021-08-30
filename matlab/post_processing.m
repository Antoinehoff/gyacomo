%% Retrieving max polynomial degree and sampling info
Npe = numel(Pe); Nje = numel(Je); [JE,PE] = meshgrid(Je,Pe);
Npi = numel(Pi); Nji = numel(Ji); [JI,PI] = meshgrid(Ji,Pi);
Ns5D      = numel(Ts5D);
Ns3D      = numel(Ts3D);
% renaming and reshaping quantity of interest
Ts5D      = Ts5D';
Ts3D      = Ts3D';

%% Build grids
Nkx = numel(kx); Nky = numel(ky);
[KY,KX] = meshgrid(ky,kx);
Lkx = max(kx)-min(kx); Lky = max(ky)-min(ky);
dkx = Lkx/(Nkx-1); dky = Lky/(Nky-1);
KPERP2 = KY.^2+KX.^2;
[~,ikx0] = min(abs(kx)); [~,iky0] = min(abs(ky));
[KY_XY,KX_XY] = meshgrid(ky,kx);
[KZ_XZ,KX_XZ] = meshgrid(z,kx);
[KZ_YZ,KY_YZ] = meshgrid(z,ky);

Lk = max(Lkx,Lky);
Nx = max(Nkx,Nky); Ny = Nx;      Nz = numel(z);
dx = 2*pi/Lk;      dy = 2*pi/Lk; dz = 2*pi/Nz;
x = dx*(-Nx/2:(Nx/2-1)); Lx = max(x)-min(x);
y = dy*(-Ny/2:(Ny/2-1)); Ly = max(y)-min(y);
z = dz * (1:Nz);
[Y_XY,X_XY] = meshgrid(y,x);
[Z_XZ,X_XZ] = meshgrid(z,x);
[Z_YZ,Y_YZ] = meshgrid(z,y);

%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Analysis :')
disp('- iFFT')
% IFFT (Lower case = real space, upper case = frequency space)
ne00   = zeros(Nx,Ny,Nz,Ns3D); % Gyrocenter density
ni00   = zeros(Nx,Ny,Nz,Ns3D); % "
phi    = zeros(Nx,Ny,Nz,Ns3D); % Electrostatic potential
Z_phi  = zeros(Nx,Ny,Nz,Ns3D); % Zonal "
dens_e = zeros(Nx,Ny,Nz,Ns3D); % Particle density
dens_i = zeros(Nx,Ny,Nz,Ns3D); %"
Z_n_e  = zeros(Nx,Ny,Nz,Ns3D); % Zonal "
Z_n_i  = zeros(Nx,Ny,Nz,Ns3D); %"
temp_e = zeros(Nx,Ny,Nz,Ns3D); % Temperature
temp_i = zeros(Nx,Ny,Nz,Ns3D); % "
Z_T_e  = zeros(Nx,Ny,Nz,Ns3D); % Zonal "
Z_T_i  = zeros(Nx,Ny,Nz,Ns3D); %"
dyTe   = zeros(Nx,Ny,Nz,Ns3D); % Various derivatives
dyTi   = zeros(Nx,Ny,Nz,Ns3D); % "
dyni   = zeros(Nx,Ny,Nz,Ns3D); % "
dxphi  = zeros(Nx,Ny,Nz,Ns3D); % "
dyphi  = zeros(Nx,Ny,Nz,Ns3D); % "
dx2phi = zeros(Nx,Ny,Nz,Ns3D); % "

for it = 1:numel(Ts3D)
    for iz = 1:numel(z)
        NE_ = Ne00(:,:,iz,it); NI_ = Ni00(:,:,iz,it); PH_ = PHI(:,:,iz,it);
        ne00  (:,:,iz,it) = real(fftshift(ifft2((NE_),Nx,Ny)));
        ni00  (:,:,iz,it) = real(fftshift(ifft2((NI_),Nx,Ny)));
        phi   (:,:,iz,it) = real(fftshift(ifft2((PH_),Nx,Ny)));
        Z_phi (:,:,iz,it) = real(fftshift(ifft2((PH_.*(KY==0)),Nx,Ny)));
        dxphi (:,:,iz,it) = real(fftshift(ifft2(1i*KX.*(PH_),Nx,Ny)));
        dx2phi(:,:,iz,it) = real(fftshift(ifft2(-KX.^2.*(PH_),Nx,Ny)));
        dyphi (:,:,iz,it) = real(fftshift(ifft2(1i*KY.*(PH_),Nx,Ny)));
        if(W_DENS)
        DENS_E_ = DENS_E(:,:,iz,it); DENS_I_ = DENS_I(:,:,iz,it);
        dyni   (:,:,iz,it) = real(fftshift(ifft2(1i*KY.*(DENS_I_),Nx,Ny)));
        dens_e (:,:,iz,it) = real(fftshift(ifft2((DENS_E_),Nx,Ny)));
        dens_i (:,:,iz,it) = real(fftshift(ifft2((DENS_I_),Nx,Ny)));
        Z_n_e  (:,:,iz,it) = real(fftshift(ifft2((DENS_E_.*(KY==0)),Nx,Ny)));
        Z_n_i  (:,:,iz,it) = real(fftshift(ifft2((DENS_I_.*(KY==0)),Nx,Ny)));
        end
        if(W_TEMP)
        TEMP_E_ = TEMP_E(:,:,iz,it); TEMP_I_ = TEMP_I(:,:,iz,it);
        dyTe(:,:,iz,it)  = real(fftshift(ifft2(1i*KY.*(TEMP_E_),Nx,Ny)));
        dyTi(:,:,iz,it)  = real(fftshift(ifft2(1i*KY.*(TEMP_I_),Nx,Ny)));
        temp_e (:,:,iz,it) = real(fftshift(ifft2((TEMP_E_),Nx,Ny)));
        temp_i (:,:,iz,it) = real(fftshift(ifft2((TEMP_I_),Nx,Ny)));
        Z_T_e  (:,:,iz,it) = real(fftshift(ifft2((TEMP_E_.*(KY==0)),Nx,Ny)));
        Z_T_i  (:,:,iz,it) = real(fftshift(ifft2((TEMP_I_.*(KY==0)),Nx,Ny)));
        end      
    end
end


% Post processing
disp('- post processing')
% particle flux
Gamma_x= zeros(Nx,Ny,Nz,Ns3D); % Radial particle transport

phi_maxx_maxy  = zeros(Nz,Ns3D);        % Time evol. of the norm of phi
phi_avgx_maxy  = zeros(Nz,Ns3D);        % Time evol. of the norm of phi
phi_maxx_avgy  = zeros(Nz,Ns3D);        % Time evol. of the norm of phi
phi_avgx_avgy  = zeros(Nz,Ns3D);        % Time evol. of the norm of phi

shear_maxx_maxy  = zeros(Nz,Ns3D);    % Time evol. of the norm of shear
shear_avgx_maxy  = zeros(Nz,Ns3D);    % Time evol. of the norm of shear
shear_maxx_avgy  = zeros(Nz,Ns3D);    % Time evol. of the norm of shear
shear_avgx_avgy  = zeros(Nz,Ns3D);    % Time evol. of the norm of shear

Ne_norm  = zeros(Pe_max,Je_max,Ns5D);  % Time evol. of the norm of Napj
Ni_norm  = zeros(Pi_max,Ji_max,Ns5D);  % .

% Kperp spectrum interpolation
%full kperp points
kperp  = reshape(sqrt(KX.^2+KY.^2),[numel(KX),1]);
% interpolated kperps
nk_noAA = floor(2/3*numel(kx));
kp_ip = kx;
[thg, rg] = meshgrid(linspace(0,pi,2*nk_noAA),kp_ip);
[xn,yn] = pol2cart(thg,rg);
[ky_s, sortIdx] = sort(ky);
[xc,yc] = meshgrid(ky_s,kx);
phi_kp_t = zeros(numel(kp_ip),Nz,Ns3D);
%
for it = 1:numel(Ts3D) % Loop over 2D aX_XYays
    for iz = 1:numel(z)
    NE_ = Ne00(:,:,iz,it); NI_ = Ni00(:,:,iz,it); PH_ = PHI(:,:,iz,it);
    phi_maxx_maxy(iz,it)   =  max( max(squeeze(phi(:,:,iz,it))));
    phi_avgx_maxy(iz,it)   =  max(mean(squeeze(phi(:,:,iz,it))));
    phi_maxx_avgy(iz,it)   = mean( max(squeeze(phi(:,:,iz,it))));
    phi_avgx_avgy(iz,it)   = mean(mean(squeeze(phi(:,:,iz,it))));
    
    if(W_DENS)
    Gamma_x(:,:,iz,it) = dens_i(:,:,iz,it).*dyphi(:,:,iz,it);
    end

    shear_maxx_maxy(iz,it)  =  max( max(squeeze(-(dx2phi(:,:,iz,it)))));
    shear_avgx_maxy(iz,it)  =  max(mean(squeeze(-(dx2phi(:,:,iz,it)))));
    shear_maxx_avgy(iz,it)  = mean( max(squeeze(-(dx2phi(:,:,iz,it)))));
    shear_avgx_avgy(iz,it)  = mean(mean(squeeze(-(dx2phi(:,:,iz,it)))));

    Z_rth = interp2(xc,yc,squeeze(mean((abs(PHI(:,sortIdx,iz,it))).^2,3)),xn,yn);
    phi_kp_t(:,iz,it) = mean(Z_rth,2);
    end
end
%
for it = 1:numel(Ts5D) % Loop over 5D aX_XYays
    [~, it2D] = min(abs(Ts3D-Ts5D(it)));
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4)/Nkx/Nky;
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4)/Nkx/Nky;
end

%% Compute primary instability growth rate
disp('- growth rate')
% Find max value of transport (end of linear mode)
[tmp,tmax] = max(GGAMMA_RI*(2*pi/Nx/Ny)^2);
[~,itmax]  = min(abs(Ts3D-tmax));
tstart = 0.1 * Ts3D(itmax); tend = 0.5 * Ts3D(itmax);
[~,its3D_lin] = min(abs(Ts3D-tstart));
[~,ite3D_lin]   = min(abs(Ts3D-tend));

g_I          = zeros(Nkx,Nky,Nz);
for ikx = 1:Nkx
    for iky = 1:Nky
        for iz = 1:Nz
            [g_I(ikx,iky,iz), ~] = LinearFit_s(Ts3D(its3D_lin:ite3D_lin),squeeze(abs(Ni00(ikx,iky,iz,its3D_lin:ite3D_lin))));
        end
    end
end
[gmax_I,ikmax_I] = max(max(g_I(1,:,:),[],2),[],3);
kmax_I = abs(ky(ikmax_I));
Bohm_transport = ETAB/ETAN*gmax_I/kmax_I^2;

%% Compute secondary instability growth rate
disp('- growth rate')
% Find max value of transport (end of linear mode)
% [tmp,tmax] = max(GGAMMA_RI*(2*pi/Nx/Ny)^2);
% [~,itmax]  = min(abs(Ts2D-tmax));
% tstart = Ts2D(itmax); tend = 1.5*Ts2D(itmax);
[~,its3D_lin] = min(abs(Ts3D-tstart));
[~,ite3D_lin]   = min(abs(Ts3D-tend));

g_II          = zeros(Nkx,Nky);
for ikx = 1:Nkx
    for iky = 1
        for iz = 1:Nz
            [g_II(ikx,iky,iz), ~] = LinearFit_s(Ts3D(its3D_lin:ite3D_lin),squeeze(abs(Ni00(ikx,iky,iz,its3D_lin:ite3D_lin))));
        end
    end
end
[gmax_II,ikmax_II] = max(max(g_II(1,:,:),[],2),[],3);
kmax_II = abs(kx(ikmax_II));
