addpath(genpath('../matlab')) % ... add
for i_ = 1
% for ETA_ =[0.6:0.1:0.9]
%% Load results
if 1% Local results
outfile ='';
outfile ='HD_study/150x75_L_100_P_2_J_1_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_3e-02';
% outfile ='HD_study/100x50_L_50_P_2_J_1_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_1e-02';
% outfile ='kobayashi/100x50_L_50_P_2_J_1_eta_0.71429_nu_1e-02_PAGK_CLOS_0_mu_0e+00';
% outfile ='v2.7_P_2_J_1/100x50_L_200_P_2_J_1_eta_0.6_nu_1e+00_SGGK_CLOS_0_mu_0e+00';
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
outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e+00_SGGK_CLOS_0_mu_1e-02/out.txt';
% outfile = outcl{i_};
% load_marconi(outfile);
BASIC.RESDIR      = ['../',outfile(46:end-8),'/'];
BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/',outfile(46:end-8),'/'];
end

%%
% JOBNUM = 1; load_results;
JOBNUMMAX = 20; compile_results %Compile the results from first output found to JOBNUMMAX if existing

%% Retrieving max polynomial degree and sampling info
Npe = numel(Pe); Nje = numel(Je); [JE,PE] = meshgrid(Je,Pe);
Npi = numel(Pi); Nji = numel(Ji); [JI,PI] = meshgrid(Ji,Pi);
Ns5D      = numel(Ts5D);
Ns2D      = numel(Ts2D);
% renaming and reshaping quantity of interest
Ts5D      = Ts5D';
Ts2D      = Ts2D';

%% Build grids
Nkr = numel(kr); Nkz = numel(kz);
[KZ,KR] = meshgrid(kz,kr);
Lkr = max(kr)-min(kr); Lkz = max(kz)-min(kz);
dkr = Lkr/(Nkr-1); dkz = Lkz/(Nkz-1);
KPERP2 = KZ.^2+KR.^2;
[~,ikr0] = min(abs(kr)); [~,ikz0] = min(abs(kz));

Lk = max(Lkr,Lkz);
dr = 2*pi/Lk; dz = 2*pi/Lk;
Nr = max(Nkr,Nkz);         Nz = Nr;
r = dr*(-Nr/2:(Nr/2-1)); Lr = max(r)-min(r);
z = dz*(-Nz/2:(Nz/2-1)); Lz = max(z)-min(z);
[ZZ,RR] = meshgrid(z,r);

%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Analysis :')
disp('- iFFT')
% IFFT (Lower case = real space, upper case = frequency space)
ne00   = zeros(Nr,Nz,Ns2D);         % Gyrocenter density
ni00   = zeros(Nr,Nz,Ns2D);
dzTe   = zeros(Nr,Nz,Ns2D);
dzTi   = zeros(Nr,Nz,Ns2D);
dzni   = zeros(Nr,Nz,Ns2D);
np_i   = zeros(Nr,Nz,Ns5D); % Ion particle density
si00   = zeros(Nr,Nz,Ns5D);
phi    = zeros(Nr,Nz,Ns2D);
dens_e = zeros(Nr,Nz,Ns2D);
dens_i = zeros(Nr,Nz,Ns2D);
temp_e = zeros(Nr,Nz,Ns2D);
temp_i = zeros(Nr,Nz,Ns2D);
drphi  = zeros(Nr,Nz,Ns2D);
dzphi  = zeros(Nr,Nz,Ns2D);
dr2phi = zeros(Nr,Nz,Ns2D);
E_turb = zeros(1,Ns2D);    % Time evol. of the turbulence energy (Pred in Kobayashi 2015)
E_ZF   = zeros(1,Ns2D);    % Time evol. of the ZF energy (Pray in Kobayashi 2015)
for it = 1:numel(Ts2D)

end
for it = 1:numel(Ts2D)
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    ne00(:,:,it)    = real(fftshift(ifft2((NE_),Nr,Nz)));
    ni00(:,:,it)    = real(fftshift(ifft2((NI_),Nr,Nz)));
    phi (:,:,it)    = real(fftshift(ifft2((PH_),Nr,Nz)));
    drphi(:,:,it) = real(fftshift(ifft2(1i*KR.*(PH_),Nr,Nz)));
    dr2phi(:,:,it)= real(fftshift(ifft2(-KR.^2.*(PH_),Nr,Nz)));
    dzphi(:,:,it) = real(fftshift(ifft2(1i*KZ.*(PH_),Nr,Nz)));
    E_turb(it) = sum(sum((1+KR.^2+KZ.^2).*abs(PHI(:,:,it)).^2))- sum((1+kr.^2).*abs(PHI(:,1,it)).^2);
    E_ZF(it)   = kr(ikZF)^2*abs(PHI(ikZF,1,it)).^2;
    if(W_DENS && W_TEMP)
    DENS_E_ = DENS_E(:,:,it); DENS_I_ = DENS_I(:,:,it);
    TEMP_E_ = TEMP_E(:,:,it); TEMP_I_ = TEMP_I(:,:,it);
    dzni(:,:,it)  = real(fftshift(ifft2(1i*KZ.*(DENS_I_),Nr,Nz)));
    dzTe(:,:,it)  = real(fftshift(ifft2(1i*KZ.*(TEMP_E_),Nr,Nz)));
    dzTi(:,:,it)  = real(fftshift(ifft2(1i*KZ.*(TEMP_I_),Nr,Nz)));
    dens_e (:,:,it) = real(fftshift(ifft2((DENS_E_),Nr,Nz)));
    dens_i (:,:,it) = real(fftshift(ifft2((DENS_I_),Nr,Nz)));
    temp_e (:,:,it) = real(fftshift(ifft2((TEMP_E_),Nr,Nz)));
    temp_i (:,:,it) = real(fftshift(ifft2((TEMP_I_),Nr,Nz)));
    end
end

% Building a version of phi only 5D sampling times
PHI_Ts5D = zeros(Nkr,Nkz,Ns5D);
err = 0;
for it = 1:numel(Ts5D) % Loop over 5D arrays
    [shift, it2D] = min(abs(Ts2D-Ts5D(it)));
    if shift > err; err = shift; end;
    PHI_Ts5D(:,:,it) = PHI(:,:,it2D);
end
if err > 0; disp('WARNING Ts2D and Ts5D are shifted'); end;

Np_i = zeros(Nkr,Nkz,Ns5D); % Ion particle density in Fourier space

for it = 1:numel(Ts5D)
    [~, it2D] = min(abs(Ts2D-Ts5D(it)));    
    Np_i(:,:,it) = 0;
    for ij = 1:Nji
        Kn = (KPERP2/2.).^(ij-1) .* exp(-KPERP2/2)/(factorial(ij-1));
        Np_i(:,:,it) = Np_i(:,:,it) + Kn.*squeeze(Nipj(1,ij,:,:,it));
    end
    
    np_i(:,:,it)      = real(fftshift(ifft2(squeeze(Np_i(:,:,it)),Nr,Nz)));
end

% Post processing
disp('- post processing')
% gyrocenter and particle flux from real space
GFlux_ri  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ni drphi>
GFlux_zi  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ni dzphi>
GFlux_re  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ne drphi>
GFlux_ze  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ne dzphi>
PFlux_ri  = zeros(1,Ns5D);      % Particle   flux
% gyrocenter and particle flux from fourier coefficients
GFLUX_RI = real(squeeze(sum(sum(-1i*KZ.*Ni00.*conj(PHI),1),2)))*(2*pi/Nr/Nz)^2;
PFLUX_RI = real(squeeze(sum(sum(-1i*KZ.*Np_i.*conj(PHI_Ts5D),1),2)))*(2*pi/Nr/Nz)^2;
% Heat flux
Q_RI = -squeeze(mean(mean(dzphi.*temp_i,1),2))';
% Hermite energy spectrum
epsilon_e_pj = zeros(Pe_max,Je_max,Ns5D);
epsilon_i_pj = zeros(Pi_max,Ji_max,Ns5D);

phi_maxr_maxz  = zeros(1,Ns2D);        % Time evol. of the norm of phi
phi_avgr_maxz  = zeros(1,Ns2D);        % Time evol. of the norm of phi
phi_maxr_avgz  = zeros(1,Ns2D);        % Time evol. of the norm of phi
phi_avgr_avgz  = zeros(1,Ns2D);        % Time evol. of the norm of phi

shear_maxr_maxz  = zeros(1,Ns2D);    % Time evol. of the norm of shear
shear_avgr_maxz  = zeros(1,Ns2D);    % Time evol. of the norm of shear
shear_maxr_avgz  = zeros(1,Ns2D);    % Time evol. of the norm of shear
shear_avgr_avgz  = zeros(1,Ns2D);    % Time evol. of the norm of shear

Ne_norm  = zeros(Pe_max,Je_max,Ns5D);  % Time evol. of the norm of Napj
Ni_norm  = zeros(Pi_max,Ji_max,Ns5D);  % .

Ddr = 1i*KR; Ddz = 1i*KZ; lapl   = Ddr.^2 + Ddz.^2; 
% Kperp spectrum interpolation
%full kperp points
kperp  = reshape(sqrt(KR.^2+KZ.^2),[numel(KR),1]);
% interpolated kperps
nk_noAA = floor(2/3*numel(kr));
kp_ip = kr;
[thg, rg] = meshgrid(linspace(0,pi,2*nk_noAA),kp_ip);
[xn,yn] = pol2cart(thg,rg);
[kz_s, sortIdx] = sort(kz);
[xc,yc] = meshgrid(kz_s,kr);
phi_kp_t = zeros(numel(kp_ip),Ns2D);
%
for it = 1:numel(Ts2D) % Loop over 2D arrays
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    phi_maxr_maxz(it)   =  max( max(squeeze(phi(:,:,it))));
    phi_avgr_maxz(it)   =  max(mean(squeeze(phi(:,:,it))));
    phi_maxr_avgz(it)   = mean( max(squeeze(phi(:,:,it))));
    phi_avgr_avgz(it)   = mean(mean(squeeze(phi(:,:,it))));

    shear_maxr_maxz(it)  =  max( max(squeeze(-(dr2phi(:,:,it)))));
    shear_avgr_maxz(it)  =  max(mean(squeeze(-(dr2phi(:,:,it)))));
    shear_maxr_avgz(it)  = mean( max(squeeze(-(dr2phi(:,:,it)))));
    shear_avgr_avgz(it)  = mean(mean(squeeze(-(dr2phi(:,:,it)))));

    GFlux_ri(it)  = sum(sum(ni00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_zi(it)  = sum(sum(-ni00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_re(it)  = sum(sum(ne00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_ze(it)  = sum(sum(-ne00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
    
    Z_rth = interp2(xc,yc,squeeze(mean((abs(PHI(:,sortIdx,it))).^2,3)),xn,yn);
    phi_kp_t(:,it) = mean(Z_rth,2);
end
%
for it = 1:numel(Ts5D) % Loop over 5D arrays
    [~, it2D] = min(abs(Ts2D-Ts5D(it)));
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    epsilon_e_pj(:,:,it) = sqrt(pi)/2*sum(sum(abs(Nepj(:,:,:,:,it)).^2,3),4);
    epsilon_i_pj(:,:,it) = sqrt(pi)/2*sum(sum(abs(Nipj(:,:,:,:,it)).^2,3),4);
    % Particle flux
    PFlux_ri(it)   = sum(sum(np_i(:,:,it).*dzphi(:,:,it2D)))*dr*dz/Lr/Lz;
end

%% Compute primary instability growth rate
disp('- growth rate')
% Find max value of transport (end of linear mode)
[tmp,tmax] = max(GGAMMA_RI*(2*pi/Nr/Nz)^2);
[~,itmax]  = min(abs(Ts2D-tmax));
tstart = 0.1 * Ts2D(itmax); tend = 0.5 * Ts2D(itmax);
[~,its2D_lin] = min(abs(Ts2D-tstart));
[~,ite2D_lin]   = min(abs(Ts2D-tend));

g_I          = zeros(Nkr,Nkz);
for ikr = 1:Nkr
    for ikz = 1:Nkz
        [g_I(ikr,ikz), ~] = LinearFit_s(Ts2D(its2D_lin:ite2D_lin),squeeze(abs(Ni00(ikr,ikz,its2D_lin:ite2D_lin))));
    end
end
[gmax_I,ikmax_I] = max(g_I(1,:));
kmax_I = abs(kz(ikmax_I));
Bohm_transport = ETAB/ETAN*gmax_I/kmax_I^2;

%% Compute secondary instability growth rate
disp('- growth rate')
% Find max value of transport (end of linear mode)
% [tmp,tmax] = max(GGAMMA_RI*(2*pi/Nr/Nz)^2);
% [~,itmax]  = min(abs(Ts2D-tmax));
% tstart = Ts2D(itmax); tend = 1.5*Ts2D(itmax);
[~,its2D_lin] = min(abs(Ts2D-tstart));
[~,ite2D_lin]   = min(abs(Ts2D-tend));

g_II          = zeros(Nkr,Nkz);
for ikr = 1:Nkr
    for ikz = 1
        [g_II(ikr,ikz), ~] = LinearFit_s(Ts2D(its2D_lin:ite2D_lin),squeeze(abs(Ni00(ikr,ikz,its2D_lin:ite2D_lin))));
    end
end
[gmax_II,ikmax_II] = max(g_II(1,:));
kmax_II = abs(kr(ikmax_II));

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 1
%% Time evolutions and growth rate
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM),'_',PARAMS];
set(gcf, 'Position',  [100, 100, 900, 800])
subplot(111); 
    suptitle(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
        ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
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
        plot(Ts0D,GGAMMA_RI*(2*pi/Nr/Nz)^2); hold on;
%         plot(Ts2D,GFLUX_RI)
        plot(Ts0D,PGAMMA_RI*(2*pi/Nr/Nz)^2);
%         plot(Ts5D,PFLUX_RI,'--');
        legend(['Gyro. flux';'Part. flux']);
        grid on; xlabel('$t c_s/R$'); ylabel('$\Gamma_{r,i}$')
%         ylim([0,2.0]);
    if(0)
    subplot(223)
        plot(kz,g_I(1,:),'-','DisplayName','Primar. instability'); hold on;
        plot(kr,g_II(:,1),'x-','DisplayName','Second. instability'); hold on;
        plot([max(kz)*2/3,max(kz)*2/3],[0,10],'--k', 'DisplayName','2/3 Orszag AA');
        grid on; xlabel('$k\rho_s$'); ylabel('$\gamma R/c_s$'); legend('show');
        ylim([0,max(g_I(1,:))]); xlim([0,max(kz)]);
        shearplot = 426; phiplot = 428;
    else
    shearplot = 223; phiplot = 224;      
    end
    subplot(shearplot)
        clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
        lstyle   = line_styles(min(ij,numel(line_styles)));
        plot(Ts2D,shear_maxr_maxz,'DisplayName','$\max_{r,z}(s)$'); hold on;
        plot(Ts2D,shear_maxr_avgz,'DisplayName','$\max_{r}\langle s \rangle_z$'); hold on;
        plot(Ts2D,shear_avgr_maxz,'DisplayName','$\max_{z}\langle s \rangle_r$'); hold on;
        plot(Ts2D,shear_avgr_avgz,'DisplayName','$\langle s \rangle_{r,z}$'); hold on;
    grid on; xlabel('$t c_s/R$'); ylabel('$shear$'); 
    subplot(phiplot)
        clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
        lstyle   = line_styles(min(ij,numel(line_styles)));
        plot(Ts2D,phi_maxr_maxz,'DisplayName','$\max_{r,z}(\phi)$'); hold on;
        plot(Ts2D,phi_maxr_avgz,'DisplayName','$\max_{r}\langle\phi\rangle_z$'); hold on;
        plot(Ts2D,phi_avgr_maxz,'DisplayName','$\max_{z}\langle\phi\rangle_r$'); hold on;
        plot(Ts2D,phi_avgr_avgz,'DisplayName','$\langle\phi\rangle_{r,z}$'); hold on;
    grid on; xlabel('$t c_s/R$'); ylabel('$E.S. pot$');
save_figure
end

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG = 1000; % Averaging time duration
%Compute steady radial transport
tend = Ts0D(end); tstart = tend - TAVG;
[~,its0D] = min(abs(Ts0D-tstart));
[~,ite0D]   = min(abs(Ts0D-tend));
SCALE = (2*pi/Nr/Nz)^2;
gamma_infty_avg = mean(PGAMMA_RI(its0D:ite0D))*SCALE;
gamma_infty_std = std (PGAMMA_RI(its0D:ite0D))*SCALE;
% Compute steady shearing rate
tend = Ts2D(end); tstart = tend - TAVG;
[~,its2D] = min(abs(Ts2D-tstart));
[~,ite2D]   = min(abs(Ts2D-tend));
shear_infty_avg = mean(shear_maxr_avgz(its2D:ite2D));
shear_infty_std = std (shear_maxr_avgz(its2D:ite2D));
Q_infty_avg     = mean(Q_RI(its2D:ite2D))*SCALE;
Q_infty_std     = std(Q_RI(its2D:ite2D))*SCALE;
% plots
fig = figure; FIGNAME = ['ZF_transport_drphi','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
    yyaxis left
        plot(Ts0D,PGAMMA_RI*SCALE,'DisplayName','$\langle n_i d\phi/dz \rangle_z$'); hold on;
        plot(Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*gamma_infty_avg, '-k',...
            'DisplayName',['$\Gamma^{\infty} = $',num2str(gamma_infty_avg),'$\pm$',num2str(gamma_infty_std)]);
        grid on; set(gca,'xticklabel',[]); ylabel('$\Gamma_r$')
        ylim([0,5*abs(gamma_infty_avg)]); xlim([0,Ts0D(end)]);
        title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
        ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
    yyaxis right     
        plot(Ts2D,shear_maxr_avgz,'DisplayName','$\max_{r}\langle s_\phi\rangle_z$'); hold on;
        plot(Ts2D(its2D:ite2D),ones(ite2D-its2D+1,1)*shear_infty_avg, '-k',...
        'DisplayName',['$s^{\infty} = $',num2str(shear_infty_avg),'$\pm$',num2str(shear_infty_std)]);
        ylim([0,shear_infty_avg*5.0]); xlim([0,Ts0D(end)]);
        grid on; ylabel('Shear amp.');set(gca,'xticklabel',[]);% legend('show');
    subplot(312)
    yyaxis left
        plot(Ts2D,SCALE*E_ZF);
        ylabel('ZF energy');
    yyaxis right     
        plot(Ts2D,SCALE*E_turb);
        ylabel('Turb. energy');  ylim([0;1.2*max(SCALE*E_ZF(floor(0.8*numel(Ts2D)):end))]);
    subplot(313)
        [TY,TX] = meshgrid(r,Ts2D);
%         pclr = pcolor(TX,TY,squeeze(mean(drphi(:,:,:),2))'); set(pclr, 'edgecolor','none'); legend('$\langle \partial_r\phi\rangle_z$') %colorbar;
        pclr = pcolor(TX,TY,squeeze(mean(dr2phi(:,:,:),2))'); set(pclr, 'edgecolor','none'); legend('Shear ($\langle \partial_r^2\phi\rangle_z$)') %colorbar; 
        caxis(1*shear_infty_avg*[-1 1]); xlabel('$t c_s/R$'), ylabel('$r/\rho_s$'); colormap(bluewhitered(256))
save_figure
end

if 1
%% Space time diagramms
tstart = 0; tend = Ts2D(end);
[~,itstart] = min(abs(Ts2D-tstart));
[~,itend]   = min(abs(Ts2D-tend));
trange = itstart:itend;
[TY,TX] = meshgrid(r,Ts2D(trange));
fig = figure; FIGNAME = ['space_time','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(211)
%         pclr = pcolor(TX,TY,squeeze(mean(dens_i(:,:,trange).*dzphi(:,:,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
        pclr = pcolor(TX,TY,squeeze(mean(ni00(:,:,trange).*dzphi(:,:,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
        shading interp
        colormap hot;
        caxis([0.0,0.05*max(max(mean(ni00(:,:,its2D:ite2D).*dzphi(:,:,its2D:ite2D),2)))]);
        caxis([0.0,0.05]); c = colorbar; c.Label.String ='\langle n_i\partial_z\phi\rangle_z';
         xticks([]); ylabel('$r/\rho_s$')
%         legend('Radial part. transport $\langle n_i\partial_z\phi\rangle_z$')
        title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
        ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
%     subplot(312)
%         pclr = pcolor(TX,TY,squeeze(mean(temp_i(:,:,trange).*dzphi(:,:,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
%         shading interp
% %         colormap(bluewhitered(256));
%          xticks([]); ylabel('$r/\rho_s$')
%         legend('Radial part. transport $\langle T_i\partial_z\phi\rangle_z$')
%         title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
%         ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
%         ' $\mu_{hd}=$',num2str(MU)]);
    subplot(212)
        pclr = pcolor(TX,TY,squeeze(mean(drphi(:,:,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
        fieldmax = max(max(mean(abs(drphi(:,:,its2D:ite2D)),2)));
        caxis([-fieldmax,fieldmax]); c = colorbar; c.Label.String ='\langle \partial_r\phi\rangle_z';
        xlabel('$t c_s/R$'), ylabel('$r/\rho_s$')
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
xlim([-L/2,L/2]); xlabel('$r/\rho_s$'); grid on; legend('show')
end

if 1
%% |phi_k|^2 spectra (Kobayashi 2015 fig 3)
tstart = 4000; tend = 4000;
[~,itstart] = min(abs(Ts2D-tstart));
[~,itend]   = min(abs(Ts2D-tend));
trange = itstart:itend;
%full kperp points
phi_k_2 = reshape(mean((abs(PHI(:,:,trange))).^2,3),[numel(KR),1]);
kperp  = reshape(sqrt(KR.^2+KZ.^2),[numel(KR),1]);
% interpolated kperps
nk_noAA = floor(2/3*numel(kr));
kp_ip = kr;
[thg, rg] = meshgrid(linspace(0,pi,2*nk_noAA),kp_ip);
[xn,yn] = pol2cart(thg,rg);
[kz_s, sortIdx] = sort(kz);
[xc,yc] = meshgrid(kz_s,kr);
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
title({['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
', $L=',num2str(L),'$, $N=',num2str(Nr),'$'];[' $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
' $\mu_{hd}=$',num2str(MU)]});
xlim([0.1,10]);
save_figure
clear Z_rth phi_kp ni_kp Ti_kp
end

%%
t0    =000;
[~, it02D] = min(abs(Ts2D-t0));
[~, it05D] = min(abs(Ts5D-t0));
skip_ = 1; 
DELAY = 1e-3*skip_;
FRAMES_2D = it02D:skip_:numel(Ts2D);
FRAMES_5D = it05D:skip_:numel(Ts5D);
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%% part density electron
GIFNAME = ['ne',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
FIELD = real(dens_e); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$n_e$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
% create_mov
end
if 0
%% part temperature electron
GIFNAME = ['Te',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
FIELD = real(temp_e); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$T_e$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
% create_mov
end
if 0
%% part density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
FIELD = real(dens_i); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$n_i$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
% create_mov
end
if 0
%% part temperature ion
GIFNAME = ['Ti',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
FIELD = real(temp_i); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$T_i$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
% create_mov
end
if 0
%% GC Density ion
GIFNAME = ['ni00',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
FIELD = real(ni00); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$n_i^{00}$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
% create_gif
create_mov
end
if 0
%% GC Density electrons
GIFNAME = ['ne00',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
FIELD = real(ne00); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$n_e^{00}$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
% create_gif
create_mov
end
if 0
%% Phi real space
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM),'_',PARAMS];INTERP = 0;
FIELD = real(phi); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$\phi$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
% create_mov
end
if 0
%% shear
GIFNAME = ['shear_r',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 1;
FIELD = -dr2phi; X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$s$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% phi averaged on z
GIFNAME = ['phi_z0',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0; SCALING=1;
FIELD =(squeeze(mean(real(phi),2))); linestyle = '-.k'; FRAMES = FRAMES_2D;
X = (r); T = Ts2D; YMIN = -1.1; YMAX = 1.1; XMIN = min(r); XMAX = max(r);
FIELDNAME = '$\langle\phi\rangle_{z}$'; XNAME = '$r/\rho_s$';
create_gif_1D_phi
end
if 0
%% flow averaged on z
GIFNAME = ['zf_z0',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0; SCALING=1;
FIELD =(squeeze(mean(real(drphi),2))); linestyle = '-.k'; FRAMES = FRAMES_2D;
X = (r); T = Ts2D; YMIN = -1.1; YMAX = 1.1; XMIN = min(r); XMAX = max(r);
FIELDNAME = '$\langle\phi\rangle_{z}$'; XNAME = '$r/\rho_s$';
create_gif_1D_phi
end
if 0
%% shear averaged on z
GIFNAME = ['shear_z0',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0; SCALING=1;
FIELD =(squeeze(mean(real(dr2phi),2))); linestyle = '-.k'; FRAMES = FRAMES_2D;
X = (r); T = Ts2D; YMIN = -1.1; YMAX = 1.1; XMIN = min(r); XMAX = max(r);
FIELDNAME = '$\langle\phi\rangle_{z}$'; XNAME = '$r/\rho_s$';
create_gif_1D_phi
end
if 0
%% phi kperp spectrum
GIFNAME = ['phi2_kp',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0; SCALING = 0;
FIELD =log10(phi_kp_t); linestyle = '-'; FRAMES = FRAMES_2D;
X = kp_ip; T = Ts2D; YMIN = -20; YMAX = 10; XMIN = min(kr); XMAX = max(kr);
FIELDNAME = '$\log_{10}|\tilde\phi_k|^2$'; XNAME = '$k_r\rho_s$';
create_gif_1D
end
if 0
%% Density ion frequency
GIFNAME = ['Ni00',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0; FRAMES = FRAMES_2D;
FIELD =ifftshift((abs(Ni00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts2D;
FIELDNAME = '$N_i^{00}$'; XNAME = '$k_r\rho_s$'; YNAME = '$k_z\rho_s$';
create_gif
end
if 0
%% Density electron frequency
GIFNAME = ['Ne00',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0; FRAMES = FRAMES_2D;
FIELD =ifftshift((abs(Ne00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts2D;
FIELDNAME = '$N_e^{00}$'; XNAME = '$k_r\rho_s$'; YNAME = '$k_z\rho_s$';
create_gif
end
if 0
%% kr vs P Si
GIFNAME = ['Sip0_kr',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
plt = @(x) squeeze(max((abs(x)),[],4));
FIELD =plt(Sipj(:,1,:,:,:)); X = kr'; Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = '$N_i^{p0}$'; XNAME = '$k_{max}\rho_s$'; YNAME = '$P$, $\sum_z$';
create_gif_imagesc
end
if 0
%% maxkz, kr vs p, for all Nipj over time
GIFNAME = ['Nipj_kr',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
plt = @(x) squeeze(sum((abs(x)),4));
FIELD = plt(Nipj); X = kr'; Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_i'; XNAME = '$k_r\rho_s$'; YNAME = '$P$';
create_gif_5D
end
if 0
%% maxkr, kz vs p, for all Nipj over time
GIFNAME = ['Nipj_kz',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
plt = @(x) fftshift(squeeze(sum((abs(x)),3)),3);
FIELD = plt(Nipj); X = sort(kz'); Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_i'; XNAME = '$k_z\rho_s$'; YNAME = '$P$, $\sum_r$';
create_gif_5D
end
%%
ZF_fourier_analysis
end