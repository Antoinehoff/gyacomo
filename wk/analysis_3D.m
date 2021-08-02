addpath(genpath('../matlab')) % ... add
for i_ = 1
% for ETA_ =[0.6:0.1:0.9]
%% Load results
if 1% Local results
outfile ='';
outfile ='Blob_diffusion/100x50_L_60_P_2_J_1_eta_Inf_nu_1e-01_DGGK_mu_0e+00';
% outfile ='test_3D/100x50x10_L_60_q0_1_P_2_J_1_eta_0.6_nu_1e-01_DGGK_mu_2e-03';
% outfile ='test_3D/100x50_L_60_P_2_J_1_eta_0.6_nu_1e-01_DGGK_mu_0e+00';
    BASIC.RESDIR      = ['../results/',outfile,'/'];
    BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
    CMD = ['cp ', BASIC.RESDIR,'outputs* ',BASIC.MISCDIR]; disp(CMD);
    system(CMD);
end
if 0% Marconi results
outfile ='';
outfile ='';
outfile ='';
outfile ='';
outfile ='';
outfile ='';
BASIC.RESDIR      = ['../',outfile(46:end-8),'/'];
BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/',outfile(46:end-8),'/'];
end

%%
JOBNUMMIN = 00; JOBNUMMAX = 20; 
compile_results_3D %Compile the results from first output found to JOBNUMMAX if existing

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

Lk = max(Lkx,Lky);
Nx = max(Nkx,Nky); Ny = Nx;      Nz = numel(z);
dx = 2*pi/Lk;      dy = 2*pi/Lk; dz = q0*2*pi/Nz;
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
ne00   = zeros(Nx,Ny,Nz,Ns3D);         % Gyrocenter density
ni00   = zeros(Nx,Ny,Nz,Ns3D);
dzTe   = zeros(Nx,Ny,Nz,Ns3D);
dzTi   = zeros(Nx,Ny,Nz,Ns3D);
dzni   = zeros(Nx,Ny,Nz,Ns3D);
np_i   = zeros(Nx,Ny,Nz,Ns5D); % Ion particle density
si00   = zeros(Nx,Ny,Nz,Ns5D);
phi    = zeros(Nx,Ny,Nz,Ns3D);
dens_e = zeros(Nx,Ny,Nz,Ns3D);
dens_i = zeros(Nx,Ny,Nz,Ns3D);
temp_e = zeros(Nx,Ny,Nz,Ns3D);
temp_i = zeros(Nx,Ny,Nz,Ns3D);
drphi  = zeros(Nx,Ny,Nz,Ns3D);
dzphi  = zeros(Nx,Ny,Nz,Ns3D);
dr2phi = zeros(Nx,Ny,Nz,Ns3D);

for it = 1:numel(Ts3D)
    for iz = 1:numel(z)
        NE_ = Ne00(:,:,iz,it); NI_ = Ni00(:,:,iz,it); PH_ = PHI(:,:,iz,it);
        ne00(:,:,iz,it)    = real(fftshift(ifft2((NE_),Nx,Ny)));
        ni00(:,:,iz,it)    = real(fftshift(ifft2((NI_),Nx,Ny)));
        phi (:,:,iz,it)    = real(fftshift(ifft2((PH_),Nx,Ny)));
        drphi(:,:,iz,it) = real(fftshift(ifft2(1i*KX.*(PH_),Nx,Ny)));
        dr2phi(:,:,iz,it)= real(fftshift(ifft2(-KX.^2.*(PH_),Nx,Ny)));
        dzphi(:,:,iz,it) = real(fftshift(ifft2(1i*KY.*(PH_),Nx,Ny)));
        if(W_DENS && W_TEMP)
        DENS_E_ = DENS_E(:,:,iz,it); DENS_I_ = DENS_I(:,:,iz,it);
        TEMP_E_ = TEMP_E(:,:,iz,it); TEMP_I_ = TEMP_I(:,:,iz,it);
        dzni(:,:,iz,it)  = real(fftshift(ifft2(1i*KY.*(DENS_I_),Nx,Ny)));
        dzTe(:,:,iz,it)  = real(fftshift(ifft2(1i*KY.*(TEMP_E_),Nx,Ny)));
        dzTi(:,:,iz,it)  = real(fftshift(ifft2(1i*KY.*(TEMP_I_),Nx,Ny)));
        dens_e (:,:,iz,it) = real(fftshift(ifft2((DENS_E_),Nx,Ny)));
        dens_i (:,:,iz,it) = real(fftshift(ifft2((DENS_I_),Nx,Ny)));
        temp_e (:,:,iz,it) = real(fftshift(ifft2((TEMP_E_),Nx,Ny)));
        temp_i (:,:,iz,it) = real(fftshift(ifft2((TEMP_I_),Nx,Ny)));
        end
    end
end

Np_i = zeros(Nkx,Nky,Ns5D); % Ion particle density in Fourier space

for it = 1:numel(Ts5D)
    [~, it2D] = min(abs(Ts3D-Ts5D(it)));    
    Np_i(:,:,it) = 0;
    for ij = 1:Nji
        Kn = (KPERP2/2.).^(ij-1) .* exp(-KPERP2/2)/(factorial(ij-1));
        Np_i(:,:,it) = Np_i(:,:,it) + Kn.*squeeze(Nipj(1,ij,:,:,it));
    end
    
    np_i(:,:,it)      = real(fftshift(ifft2(squeeze(Np_i(:,:,it)),Nx,Ny)));
end

% Post processing
disp('- post processing')
% gyrocenter and particle flux from real space
GFlux_xi  = zeros(1,Ns3D);      % Gyrocenter flux Gamma = <ni drphi>
GFlux_yi  = zeros(1,Ns3D);      % Gyrocenter flux Gamma = <ni dzphi>
GFlux_xe  = zeros(1,Ns3D);      % Gyrocenter flux Gamma = <ne drphi>
GFlux_ye  = zeros(1,Ns3D);      % Gyrocenter flux Gamma = <ne dzphi>
% Hermite energy spectrum
epsilon_e_pj = zeros(Pe_max,Je_max,Ns5D);
epsilon_i_pj = zeros(Pi_max,Ji_max,Ns5D);

phi_maxx_maxy  = zeros(Nz,Ns3D);        % Time evol. of the norm of phi
phi_avgx_maxy  = zeros(Nz,Ns3D);        % Time evol. of the norm of phi
phi_maxx_avg  = zeros(Nz,Ns3D);        % Time evol. of the norm of phi
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

    shear_maxx_maxy(iz,it)  =  max( max(squeeze(-(dr2phi(:,:,iz,it)))));
    shear_avgx_maxy(iz,it)  =  max(mean(squeeze(-(dr2phi(:,:,iz,it)))));
    shear_maxx_avgy(iz,it)  = mean( max(squeeze(-(dr2phi(:,:,iz,it)))));
    shear_avgx_avgy(iz,it)  = mean(mean(squeeze(-(dr2phi(:,:,iz,it)))));

    GFlux_xi(iz,it)  = sum(sum(ni00(:,:,iz,it).*dzphi(:,:,iz,it)))*dx*dy/Lx/Ly;
    GFlux_yi(iz,it)  = sum(sum(-ni00(:,:,iz,it).*drphi(:,:,iz,it)))*dx*dy/Lx/Ly;
    GFlux_xe(iz,it)  = sum(sum(ne00(:,:,iz,it).*dzphi(:,:,iz,it)))*dx*dy/Lx/Ly;
    GFlux_ye(iz,it)  = sum(sum(-ne00(:,:,iz,it).*drphi(:,:,iz,it)))*dx*dy/Lx/Ly;
    
    Z_rth = interp2(xc,yc,squeeze(mean((abs(PHI(:,sortIdx,iz,it))).^2,3)),xn,yn);
    phi_kp_t(:,iz,it) = mean(Z_rth,2);
    end
end
%
for it = 1:numel(Ts5D) % Loop over 5D aX_XYays
    [~, it2D] = min(abs(Ts3D-Ts5D(it)));
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4)/Nkx/Nky;
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4)/Nkx/Nky;
    epsilon_e_pj(:,:,it) = sqrt(pi)/2*sum(sum(abs(Nepj(:,:,:,:,it)).^2,3),4);
    epsilon_i_pj(:,:,it) = sqrt(pi)/2*sum(sum(abs(Nipj(:,:,:,:,it)).^2,3),4);
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

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 1
%% Time evolutions and growth rate
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM),'_',PARAMS];
set(gcf, 'Position',  [100, 100, 900, 800])
subplot(111); 
    suptitle(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta=$',num2str(ETAB/ETAN),...
        ', $L=',num2str(L),'$, $N=',num2str(Nx),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
    subplot(421); 
    for ip = 1:Pe_max
        for ij = 1:Je_max
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_e^{',num2str(ip-1),num2str(ij-1),'}$'];
            clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
            lstyle   = line_styles(min(ij,numel(line_styles)));
            semilogy(Ts5D,plt(Ne_norm),'DisplayName',plotname,...
                'Color',clr,'LineStyle',lstyle{1}); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_e^{pj}|$');
    subplot(423)
    for ip = 1:Pi_max
        for ij = 1:Ji_max
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
            lstyle   = line_styles(min(ij,numel(line_styles)));
            semilogy(Ts5D,plt(Ni_norm),'DisplayName',plotname,...
                'Color',clr,'LineStyle',lstyle{1}); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_i^{pj}|$'); xlabel('$t c_s/R$')
    subplot(222)
        plot(Ts0D,GGAMMA_RI*(2*pi/Nx/Ny)^2); hold on;
        plot(Ts0D,PGAMMA_RI*(2*pi/Nx/Ny)^2);
        legend(['Gyro. flux';'Part. flux']);
        grid on; xlabel('$t c_s/R$'); ylabel('$\Gamma_{r,i}$')
    if(~isnan(max(max(g_I(1,:,:)))))
    subplot(223)
        plot(ky,max(g_I(1,:,:),[],3),'-','DisplayName','Primar. instability'); hold on;
        plot(kx,max(g_II(:,1,:),[],3),'x-','DisplayName','Second. instability'); hold on;
        plot([max(ky)*2/3,max(ky)*2/3],[0,10],'--k', 'DisplayName','2/3 Orszag AA');
        grid on; xlabel('$k\rho_s$'); ylabel('$\gamma R/c_s$'); legend('show');
        ylim([0,max(max(g_I(1,:,:)))]); xlim([0,max(ky)]);
        shearplot = 426; phiplot = 428;
    else
    shearplot = 223; phiplot = 224;      
    end
    subplot(shearplot)
        plt = @(x) mean(x,1);
        clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
        lstyle   = line_styles(min(ij,numel(line_styles)));
        plot(Ts3D,plt(shear_maxx_maxy),'DisplayName','$\max_{r,z}(s)$'); hold on;
        plot(Ts3D,plt(shear_maxx_avgy),'DisplayName','$\max_{r}\langle s \rangle_z$'); hold on;
        plot(Ts3D,plt(shear_avgx_maxy),'DisplayName','$\max_{z}\langle s \rangle_r$'); hold on;
        plot(Ts3D,plt(shear_avgx_avgy),'DisplayName','$\langle s \rangle_{r,z}$'); hold on;
    grid on; xlabel('$t c_s/R$'); ylabel('$shear$'); 
    subplot(phiplot)
        clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
        lstyle   = line_styles(min(ij,numel(line_styles)));
        plot(Ts3D,plt(phi_maxx_maxy),'DisplayName','$\max_{r,z}(\phi)$'); hold on;
        plot(Ts3D,plt(phi_maxx_avg),'DisplayName','$\max_{r}\langle\phi\rangle_z$'); hold on;
        plot(Ts3D,plt(phi_avgx_maxy),'DisplayName','$\max_{z}\langle\phi\rangle_r$'); hold on;
        plot(Ts3D,plt(phi_avgx_avgy),'DisplayName','$\langle\phi\rangle_{r,z}$'); hold on;
    grid on; xlabel('$t c_s/R$'); ylabel('$E.S. pot$');
save_figure
end

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG = 5000; % Averaging time duration
%Compute steady radial transport
tend = Ts0D(end); tstart = tend - TAVG;
[~,its0D] = min(abs(Ts0D-tstart));
[~,ite0D]   = min(abs(Ts0D-tend));
SCALE = (2*pi/Nx/Ny)^2;
gamma_infty_avg = mean(PGAMMA_RI(its0D:ite0D))*SCALE;
gamma_infty_std = std (PGAMMA_RI(its0D:ite0D))*SCALE;
% Compute steady shearing rate
tend = Ts3D(end); tstart = tend - TAVG;
[~,its2D] = min(abs(Ts3D-tstart));
[~,ite2D]   = min(abs(Ts3D-tend));
shear_infty_avg = mean(mean(shear_maxx_avgy(:,its2D:ite2D),1));
shear_infty_std = std (mean(shear_maxx_avgy(:,its2D:ite2D),1));
% plots
fig = figure; FIGNAME = ['ZF_transport_drphi','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
%     yyaxis left
        plot(Ts0D,PGAMMA_RI*SCALE,'DisplayName','$\langle n_i d\phi/dz \rangle_z$'); hold on;
        plot(Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*gamma_infty_avg, '-k',...
            'DisplayName',['$\Gamma^{\infty} = $',num2str(gamma_infty_avg),'$\pm$',num2str(gamma_infty_std)]);
        grid on; set(gca,'xticklabel',[]); ylabel('$\Gamma_r$')
        ylim([0,5*abs(gamma_infty_avg)]); xlim([Ts0D(1),Ts0D(end)]);
        title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta=$',num2str(ETAB/ETAN),...
        ', $L=',num2str(L),'$, $N=',num2str(Nx),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
        %         
    subplot(312)
        clr      = line_colors(1,:);
        lstyle   = line_styles(1);
%         plt = @(x_) mean(x_,1);
        plt = @(x_) x_(1,:);
        plot(Ts3D,plt(shear_maxx_maxy),'DisplayName','$\max_{r,z}(s_\phi)$'); hold on;
        plot(Ts3D,plt(shear_maxx_avgy),'DisplayName','$\max_{r}\langle s_\phi\rangle_z$'); hold on;
        plot(Ts3D,plt(shear_avgx_maxy),'DisplayName','$\max_{z}\langle s_\phi\rangle_r$'); hold on;
        plot(Ts3D,plt(shear_avgx_avgy),'DisplayName','$\langle s_\phi\rangle_{r,z}$'); hold on;
        plot(Ts3D(its2D:ite2D),ones(ite2D-its2D+1,1)*shear_infty_avg, '-k',...
        'DisplayName',['$s^{\infty} = $',num2str(shear_infty_avg),'$\pm$',num2str(shear_infty_std)]);
        ylim([0,shear_infty_avg*5.0]); xlim([Ts0D(1),Ts0D(end)]);
        grid on; ylabel('Shear amp.');set(gca,'xticklabel',[]);% legend('show');
    subplot(313)
        [TY,TX] = meshgrid(x,Ts3D);
        pclr = pcolor(TX,TY,squeeze((mean(dr2phi(:,:,1,:),2)))'); set(pclr, 'edgecolor','none'); legend('Shear ($\langle \partial_r^2\phi\rangle_z$)') %colorbar; 
        caxis(1*shear_infty_avg*[-1 1]); xlabel('$t c_s/R$'), ylabel('$x/\rho_s$'); colormap(bluewhitered(256))
save_figure
end

if 0
%% Space time diagramms
tstart = 0; tend = Ts3D(end);
[~,itstart] = min(abs(Ts3D-tstart));
[~,itend]   = min(abs(Ts3D-tend));
trange = itstart:itend;
[TY,TX] = meshgrid(x,Ts3D(trange));
fig = figure; FIGNAME = ['space_time','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(211)
%         pclr = pcolor(TX,TY,squeeze(mean(dens_i(:,:,trange).*dzphi(:,:,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
        pclr = pcolor(TX,TY,squeeze(mean(ni00(:,:,trange).*dzphi(:,:,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
        shading interp
        colormap hot;
        caxis([0.0,0.05*max(max(mean(ni00(:,:,its2D:ite2D).*dzphi(:,:,its2D:ite2D),2)))]);
        caxis([0.0,0.01]); c = colorbar; c.Label.String ='\langle n_i\partial_z\phi\rangle_z';
         xticks([]); ylabel('$x/\rho_s$')
%         legend('Radial part. transport $\langle n_i\partial_z\phi\rangle_z$')
        title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta=$',num2str(ETAB/ETAN),...
        ', $L=',num2str(L),'$, $N=',num2str(Nx),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
    subplot(212)
        pclr = pcolor(TX,TY,squeeze(mean(drphi(:,:,1,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
        fieldmax = max(max(mean(abs(drphi(:,:,1,its2D:ite2D)),2)));
        caxis([-fieldmax,fieldmax]); c = colorbar; c.Label.String ='\langle \partial_r\phi\rangle_z';
        xlabel('$t c_s/R$'), ylabel('$x/\rho_s$')
%         legend('Zonal flow $\langle \partial_r\phi\rangle_z$')        
save_figure
end

if 0
%% Averaged shear and Reynold stress profiles
trange = its2D:ite2D;
% trange = 100:200;
figure;
plt = @(x) squeeze(mean(mean(real(x(:,:,trange)),2),3))/max(abs(squeeze(mean(mean(real(x(:,:,trange)),2),3))));
plot(r,plt(-dr2phi),'-k','DisplayName','Zonal shear'); hold on;
plot(r,plt(-drphi.*dzphi),'-','DisplayName','$\Pi_\phi$'); hold on;
% plot(r,plt(-drphi.*dzTe),'-','DisplayName','$\Pi_{Te}$'); hold on;
plot(r,plt(-drphi.*dzTi),'-','DisplayName','$\Pi_{Ti}$'); hold on;
plot(r,plt(-drphi.*dzphi-drphi.*dzTi),'-','DisplayName','$\Pi_\phi+\Pi_{Ti}$'); hold on;
% plot(r,plt(-drphi.*dzphi-drphi.*dzTi-drphi.*dzTe),'-','DisplayName','$\Pi_\phi+\Pi_{Te}+\Pi_{Ti}$'); hold on;
xlim([-L/2,L/2]); xlabel('$x/\rho_s$'); grid on; legend('show')
end

if 0
%% |phi_k|^2 spectra (Kobayashi 2015 fig 3)
tstart = 0.8*Ts3D(end); tend = Ts3D(end);
[~,itstart] = min(abs(Ts3D-tstart));
[~,itend]   = min(abs(Ts3D-tend));
trange = itstart:itend;
%full kperp points
phi_k_2 = reshape(mean(mean(abs(PHI(:,:,:,trange)),3),4).^2,[numel(KX),1]);
kperp  = reshape(sqrt(KX.^2+KY.^2),[numel(KX),1]);
% interpolated kperps
nk_noAA = floor(2/3*numel(kx));
kp_ip = kx;
[thg, rg] = meshgrid(linspace(0,pi,2*nk_noAA),kp_ip);
[xn,yn] = pol2cart(thg,rg);
[ky_s, sortIdx] = sort(ky);
[xc,yc] = meshgrid(ky_s,kx);
Z_rth = interp2(xc,yc,squeeze(mean((abs(PHI(:,sortIdx,trange))).^2,3)),xn,yn);
phi_kp = mean(Z_rth,2);
Z_rth = interp2(xc,yc,squeeze(mean((abs(Ni00(:,sortIdx,trange))).^2,3)),xn,yn);
ni00_kp = mean(Z_rth,2);
Z_rth = interp2(xc,yc,squeeze(mean((abs(Ne00(:,sortIdx,trange))).^2,3)),xn,yn);
ne00_kp = mean(Z_rth,2);
%for theorical trends
a1 = phi_kp(2)*kp_ip(2).^(13/3);
a2 = phi_kp(2)*kp_ip(2).^(3)./(1+kp_ip(2).^2).^(-2);
fig = figure; FIGNAME = ['cascade','_',PARAMS];set(gcf, 'Position',  [100, 100, 500, 500])
% scatter(kperp,phi_k_2,'.k','MarkerEdgeAlpha',0.4,'DisplayName','$|\phi_k|^2$'); hold on; grid on;
plot(kp_ip,phi_kp,'^','DisplayName','$\langle|\phi_k|^2\rangle_{k_\perp}$'); hold on;
plot(kp_ip,ni00_kp, '^','DisplayName','$\langle|N_{i,k}^{00}|^2\rangle_{k_\perp}$'); hold on;
plot(kp_ip,ne00_kp, '^','DisplayName','$\langle|N_{e,k}^{00}|^2\rangle_{k_\perp}$'); hold on;
plot(kp_ip,a1*kp_ip.^(-13/3),'-','DisplayName','$k^{-13/3}$');
plot(kp_ip,a2/100*kp_ip.^(-3)./(1+kp_ip.^2).^2,'-','DisplayName','$k^{-3}/(1+k^2)^2$');
set(gca, 'XScale', 'log');set(gca, 'YScale', 'log'); grid on
xlabel('$k_\perp \rho_s$'); legend('show','Location','southwest')
title({['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta=$',num2str(ETAB/ETAN),...
', $L=',num2str(L),'$, $N=',num2str(Nx),'$'];[' $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
' $\mu_{hd}=$',num2str(MU)]});
xlim([0.1,10]);
save_figure
clear Z_rth phi_kp ni_kp Ti_kp
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
t0    =00; iz = 1; ix = 1; iy = 1;
skip_ = 1; DELAY = 1e-2*skip_;
[~, it03D] = min(abs(Ts3D-t0)); FRAMES_3D = it03D:skip_:numel(Ts3D);
[~, it05D] = min(abs(Ts5D-t0)); FRAMES_5D = it05D:skip_:numel(Ts5D);
INTERP = 0; T = Ts3D; FRAMES = FRAMES_3D;
% Field to plot
FIELD = dens_e; NAME = 'ne';   FIELDNAME = 'n_e';
% FIELD = dens_i; NAME = 'ni';   FIELDNAME = 'n_i';
% FIELD = temp_e; NAME = 'Te';   FIELDNAME = 'n_i';
% FIELD = temp_i; NAME = 'Ti';   FIELDNAME = 'n_i';
% FIELD = ne00;   NAME = 'ne00'; FIELDNAME = 'n_e^{00}';
% FIELD = ni00;   NAME = 'ni00'; FIELDNAME = 'n_i^{00}';
% Slice
% plt = @(x) real(x(ix, :, :,:)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) real(x( :,iy, :,:)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
plt = @(x) real(x( :, :,iz,:)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

% Averaged
% plt = @(x) mean(x,1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) mean(x,2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
% plt = @(x) mean(x,3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

FIELD = squeeze(plt(FIELD));

% Naming
GIFNAME   = [NAME,sprintf('_%.2d',JOBNUM),'_',PARAMS];

% Create movie (gif or mp4)
create_gif
% create_mov
end

if 0
%% Photomaton : real space

% Chose the field to plot
% FIELD = ni00;   FNAME = 'ni00'; FIELDLTX = 'n_i^{00}';
% FIELD = ne00;   FNAME = 'ne00'; FIELDLTX = 'n_e^{00}'
% FIELD = dens_i; FNAME = 'ni';   FIELDLTX = 'n_i';
% FIELD = dens_e; FNAME = 'ne';   FIELDLTX = 'n_e';
% FIELD = temp_i; FNAME = 'Ti';   FIELDLTX = 'T_i';
% FIELD = temp_e; FNAME = 'Te';   FIELDLTX = 'T_e';
FIELD = phi; FNAME = 'phi'; FIELDLTX = '\phi';

% Chose when to plot it
tf = [0 1 2 3];

% Slice
ix = 1; iy = 1; iz = 1;
% plt = @(x,it) real(x(ix, :, :,it)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(x(ix))),')']
% plt = @(x,it) real(x( :,iy, :,it)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(y=',num2str(round(y(iy))),')']
% plt = @(x,it) real(x( :, :,iz,it)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(z(iz))),')'] 

% Averaged
% plt = @(x,it) mean(x(:,:,:,it),1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_x']
% plt = @(x,it) mean(x(:,:,:,it),2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_y']
plt = @(x,it) mean(x(:,:,:,it),3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_z'] 


%
TNAME = [];
fig = figure; FIGNAME = [FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 350])
plt_2 = @(x) x./max(max(x));
    for i_ = 1:numel(tf)
    [~,it] = min(abs(Ts3D-tf(i_))); TNAME = [TNAME,'_',num2str(Ts3D(it))];
    subplot(1,numel(tf),i_)
        DATA = plt_2(squeeze(plt(FIELD,it)));
        pclr = pcolor((X),(Y),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); caxis([-1,1]);
        xlabel(latexize(XNAME)); ylabel(latexize(YNAME));set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts3D(it)));
    end
    legend(latexize(FIELDLTX));
save_figure
end

%%
% ZF_fourier_analysis
end