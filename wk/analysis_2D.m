%% Load results
outfile ='';
if 0
    %% Load from Marconi
outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HeLaZ_v2.4_eta_0.6_nu_1e-01/200x100_L_120_P_10_J_5_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-02/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HeLaZ_v2.5_eta_0.6_nu_1e-01/200x100_L_80_P_10_J_5_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_4e-03/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HeLaZ_v2.5_eta_0.6_nu_1e+00/200x100_L_120_P_12_J_6_eta_0.6_nu_1e+00_DGGK_CLOS_0_mu_2e-02/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HeLaZ_v2.5_eta_0.6_nu_1e-01/200x100_L_120_P_12_J_6_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-02/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HeLaZ_v2.5_eta_0.6_nu_1e-01/200x100_L_120_P_12_J_6_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_1e-02/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HeLaZ_v2.5_eta_0.6_nu_1e-01/200x100_L_120_P_12_J_6_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-02/out.txt';
    BASIC.RESDIR = load_marconi(outfile);
end
if 0
    %% Load from Daint
    outfile ='/scratch/snx3000/ahoffman/HeLaZ/results/HeLaZ_v2.5_eta_0.6_nu_1e-01/200x100_L_120_P_12_J_6_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-02/out.txt';
    BASIC.RESDIR = load_daint(outfile);
end
%%
% JOBNUM = 1; load_results;
compile_results

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
np_i   = zeros(Nr,Nz,Ns5D); % Ion particle density
si00   = zeros(Nr,Nz,Ns5D);
phi    = zeros(Nr,Nz,Ns2D);
drphi  = zeros(Nr,Nz,Ns2D);
dr2phi = zeros(Nr,Nz,Ns2D);
dzphi  = zeros(Nr,Nz,Ns2D);

for it = 1:numel(Ts2D)
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    ne00(:,:,it)  = real(fftshift(ifft2((NE_),Nr,Nz)));
    ni00(:,:,it)  = real(fftshift(ifft2((NI_),Nr,Nz)));
    phi (:,:,it)  = real(fftshift(ifft2((PH_),Nr,Nz)));
    drphi(:,:,it) = real(fftshift(ifft2(1i*KR.*(PH_),Nr,Nz)));
    dr2phi(:,:,it)= real(fftshift(ifft2(-KR.^2.*(PH_),Nr,Nz)));
    dzphi(:,:,it) = real(fftshift(ifft2(1i*KZ.*(PH_),Nr,Nz)));
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
E_pot    = zeros(1,Ns2D);      % Potential energy n^2
E_kin    = zeros(1,Ns2D);      % Kinetic energy grad(phi)^2
ExB      = zeros(1,Ns2D);      % ExB drift intensity \propto |\grad \phi|
% gyrocenter and particle flux from real space
GFlux_ri  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ni drphi>
GFlux_zi  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ni dzphi>
GFlux_re  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ne drphi>
GFlux_ze  = zeros(1,Ns2D);      % Gyrocenter flux Gamma = <ne dzphi>
PFlux_ri  = zeros(1,Ns5D);      % Particle   flux
% gyrocenter and particle flux from fourier coefficients
GFLUX_RI = real(squeeze(sum(sum(-1i*KZ.*Ni00.*conj(PHI),1),2)))*(2*pi/Nr/Nz)^2;
PFLUX_RI = real(squeeze(sum(sum(-1i*KZ.*Np_i.*conj(PHI_Ts5D),1),2)))*(2*pi/Nr/Nz)^2;
% Hermite energy spectrum
epsilon_e_pj = zeros(Npe,Nje,Ns5D);
epsilon_i_pj = zeros(Npi,Nji,Ns5D);

phi_maxr_maxz  = zeros(1,Ns2D);        % Time evol. of the norm of phi
phi_avgr_maxz  = zeros(1,Ns2D);        % Time evol. of the norm of phi
phi_maxr_avgz  = zeros(1,Ns2D);        % Time evol. of the norm of phi
phi_avgr_avgz  = zeros(1,Ns2D);        % Time evol. of the norm of phi
Ne_norm  = zeros(Npe,Nje,Ns5D);  % Time evol. of the norm of Napj
Ni_norm  = zeros(Npi,Nji,Ns5D);  % .

Ddr = 1i*KR; Ddz = 1i*KZ; lapl   = Ddr.^2 + Ddz.^2; 

for it = 1:numel(Ts2D) % Loop over 2D arrays
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    phi_maxr_maxz(it)   =  max( max(squeeze(phi(:,:,it))));
    phi_avgr_maxz(it)   =  max(mean(squeeze(phi(:,:,it))));
    phi_maxr_avgz(it)   = mean( max(squeeze(phi(:,:,it))));
    phi_avgr_avgz(it)   = mean(mean(squeeze(phi(:,:,it))));
    ExB(it)       = max(max(max(abs(phi(3:end,:,it)-phi(1:end-2,:,it))/(2*dr))),max(max(abs(phi(:,3:end,it)-phi(:,1:end-2,it))'/(2*dz))));
    GFlux_ri(it)  = sum(sum(ni00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_zi(it)  = sum(sum(-ni00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_re(it)  = sum(sum(ne00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_ze(it)  = sum(sum(-ne00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
end

for it = 1:numel(Ts5D) % Loop over 5D arrays
    [~, it2D] = min(abs(Ts2D-Ts5D(it)));
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    epsilon_e_pj(:,:,it) = sqrt(pi)/2*sum(sum(abs(Nepj(:,:,:,:,it)).^2,3),4);
    epsilon_i_pj(:,:,it) = sqrt(pi)/2*sum(sum(abs(Nipj(:,:,:,:,it)).^2,3),4);
    % Particle flux
    PFlux_ri(it)   = sum(sum(np_i(:,:,it).*dzphi(:,:,it2D)))*dr*dz/Lr/Lz;
end

%% Compute growth rate
disp('- growth rate')
% Find max value of transport (end of linear mode)
[tmp,tmax] = max(GGAMMA_RI*(2*pi/Nr/Nz)^2);
[~,itmax]  = min(abs(Ts2D-tmax));
tstart = 0.1 * Ts2D(itmax); tend = 0.5 * Ts2D(itmax);
g_          = zeros(Nkr,Nkz);
for ikr = 1:Nkr
    for ikz = 1:Nkz
        g_(ikr,ikz) = LinearFit_s(Ts2D,squeeze(abs(Ni00(ikr,ikz,:))),tstart,tend);
    end
end
[gmax,ikzmax] = max(g_(1,:));
kzmax = abs(kz(ikzmax));
Bohm_transport = ETAB/ETAN*gmax/kzmax^2;
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
        ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$']);
    subplot(421); 
    for ip = 1:Npe
        for ij = 1:Nje
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_e^{',num2str(Pe(ip)),num2str(Je(ij)),'}$'];
            clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
            lstyle   = line_styles(min(ij,numel(line_styles)));
            semilogy(Ts5D,plt(Ne_norm),'DisplayName',plotname,...
                'Color',clr,'LineStyle',lstyle{1}); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_e^{pj}|$');
    subplot(423)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_i^{',num2str(Pi(ip)),num2str(Ji(ij)),'}$'];
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
    subplot(223)
        plot(kz,g_(1,:),'-','DisplayName','Linear growth rate'); hold on;
        plot([max(kz)*2/3,max(kz)*2/3],[0,10],'--k', 'DisplayName','2/3 Orszag AA');
        grid on; xlabel('$k_z\rho_s$'); ylabel('$\gamma R/c_s$'); legend('show');
%         ylim([0,max(g_(1,:))]); xlim([0,max(kz)]);
    subplot(224)
        clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
        lstyle   = line_styles(min(ij,numel(line_styles)));
        plot(Ts2D,phi_maxr_maxz,'DisplayName','$\max_{r,z}(\phi)$'); hold on;
        plot(Ts2D,phi_maxr_avgz,'DisplayName','$\max_{r}\langle\phi\rangle_z$'); hold on;
        plot(Ts2D,phi_avgr_maxz,'DisplayName','$\max_{z}\langle\phi\rangle_r$'); hold on;
        plot(Ts2D,phi_avgr_avgz,'DisplayName','$\langle\phi\rangle_{r,z}$'); hold on;
    grid on; xlabel('$t c_s/R$'); ylabel('$T_e/e$'); legend('show');
save_figure
end

if 0
%% Particle fluxes
SCALING = Nkr*dkr * Nkz*dkz;
fig = figure; FIGNAME = ['gamma',sprintf('_%.2d',JOBNUM),'_',PARAMS];
set(gcf, 'Position',  [100, 100, 800, 300])
        semilogy(Ts2D,GFLUX_RI, 'color', line_colors(2,:)); hold on
        plot(Ts5D,PFLUX_RI,'.', 'color', line_colors(2,:)); hold on
        plot(Ts2D,SCALING*GFlux_ri, 'color', line_colors(1,:)); hold on
        plot(Ts5D,SCALING*PFlux_ri,'.', 'color', line_colors(1,:)); hold on
        xlabel('$tc_{s0}/\rho_s$'); ylabel('$\Gamma_r$'); grid on
        title(['$\eta=',num2str(ETAB),'\quad',...
            '\nu_{',CONAME,'}=',num2str(NU),'$', ' $P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$'])
        legend('Gyro Flux','Particle flux', 'iFFT GFlux', 'iFFT PFlux')%'$\eta\gamma_{max}/k_{max}^2$')
save_figure
end

if 0
%% Space time diagramm (fig 11 Ivanov 2020)
fig = figure; FIGNAME = ['space_time_drphi','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
        semilogy(Ts2D,GFLUX_RI,'-'); hold on
        plot(Ts5D,PFLUX_RI,'.'); hold on
%         plot(Ts2D,Bohm_transport*ones(size(Ts2D)),'--'); hold on
        ylabel('$\Gamma_r$'); grid on
        title(['$\eta=',num2str(ETAB),'\quad',...
            '\nu_{',CONAME,'}=',num2str(NU),'$'])
        legend(['$P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$'],'Particle flux')%'$\eta\gamma_{max}/k_{max}^2$')
%         set(gca,'xticklabel',[])
    subplot(312)
        yyaxis left
        plot(Ts2D,squeeze(max(max((phi)))))
        ylabel('$\max \phi$')
        yyaxis right
        plot(Ts2D,squeeze(mean(max(dr2phi))))
        ylabel('$s\sim\langle\partial_r^2\phi\rangle_z$'); grid on  
%         set(gca,'xticklabel',[])
    subplot(313)
        [TY,TX] = meshgrid(r,Ts2D);
        pclr = pcolor(TX,TY,squeeze(mean(drphi(:,:,:),2))'); set(pclr, 'edgecolor','none'); %colorbar;
        xlabel('$t c_s/R$'), ylabel('$r/\rho_s$')
        legend('$\langle\partial_r \phi\rangle_z$')
save_figure
end

if 0
%% Photomaton : real space
% FIELD = ni00; FNAME = 'ni'; XX = RR; YY = ZZ;
FIELD = phi; FNAME = 'phi'; XX = RR; YY = ZZ;
% FIELD = fftshift(abs(Ni00),2); FNAME = 'Fni'; XX = fftshift(KR,2); YY = fftshift(KZ,2);
% FIELD = fftshift(abs(PHI),2);  FNAME = 'Fphi'; XX = fftshift(KR,2); YY = fftshift(KZ,2);
tf = 100;  [~,it1] = min(abs(Ts2D-tf));
tf = 118;  [~,it2] = min(abs(Ts2D-tf)); 
tf = 140; [~,it3] = min(abs(Ts2D-tf));
tf = 300; [~,it4] = min(abs(Ts2D-tf));
fig = figure; FIGNAME = [FNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 400])
plt = @(x) x;%./max(max(x));
    subplot(141)
        DATA = plt(FIELD(:,:,it1));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$'); ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it1)));
    subplot(142)
        DATA = plt(FIELD(:,:,it2));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$'); set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it2)));
    subplot(143)
        DATA = plt(FIELD(:,:,it3));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it3)));
    subplot(144)
        DATA = plt(FIELD(:,:,it4));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$'); set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it4)));
% suptitle(['$\',FNAME,'$, $\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
%     ', $P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$']);
save_figure
end

%%
if 0
%% Show frame in kspace
tf = 0; [~,it2] = min(abs(Ts2D-tf)); [~,it5] = min(abs(Ts5D-tf));
fig = figure; FIGNAME = ['krkz_',sprintf('t=%.0f',Ts2D(it2)),'_',PARAMS];set(gcf, 'Position',  [100, 100, 700, 600])
    subplot(221); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(PHI(:,:,it2))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('$t c_s/R=%.0f$',Ts2D(it2))); legend('$|\hat\phi|$');
    subplot(222); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ni00(:,:,it2))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_i^{00}|$');
    subplot(223); plt = @(x) fftshift((abs(x)),2); FIELD = squeeze(Nipj(1,2,:,:,:));
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(FIELD(:,:,it5))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_i^{pj=01}|$');
    subplot(224); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Si00(:,:,it5))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_i^{00}$');
save_figure
end

%%
if 0
%% Hermite energy spectra
% tf = Ts2D(end-3);
fig = figure; FIGNAME = ['hermite_spectrum_',PARAMS];set(gcf, 'Position',  [100, 100, 1800, 600]);
plt = @(x) squeeze(x);
for ij = 1:Nji
    subplotnum = 100+Nji*10+ij;
    subplot(subplotnum)
    for it5 = 1:2:Ns5D
        alpha = it5*1.0/Ns5D;
        loglog(Pi(1:2:end),plt(epsilon_i_pj(1:2:end,ij,it5)),...
            'color',(1-alpha)*[0.8500, 0.3250, 0.0980]+alpha*[0, 0.4470, 0.7410],...
            'DisplayName',['t=',num2str(Ts5D(it5))]); hold on;
    end
    grid on;
    xlabel('$p$');
    TITLE = ['$\sum_{kr,kz} |N_i^{p',num2str(Ji(ij)),'}|^2$']; title(TITLE);
end
save_figure
end
%%
if 0
%% Laguerre energy spectra
% tf = Ts2D(end-3);
fig = figure; FIGNAME = ['laguerre_spectrum_',PARAMS];set(gcf, 'Position',  [100, 100, 500, 400]);
plt = @(x) squeeze(x);
for it5 = 1:2:Ns5D
    alpha = it5*1.0/Ns5D;
    loglog(Ji,plt(max(epsilon_i_pj(:,:,it5),[],1)),...
        'color',(1-alpha)*[0.8500, 0.3250, 0.0980]+alpha*[0, 0.4470, 0.7410],...
        'DisplayName',['t=',num2str(Ts5D(it5))]); hold on;
end
grid on;
xlabel('$j$');
TITLE = ['$\max_p\sum_{kr,kz} |N_i^{pj}|^2$']; title(TITLE);
save_figure
end

%%
no_AA     = (2:floor(2*Nkr/3));
tKHI      = 100;
[~,itKHI] = min(abs(Ts2D-tKHI));
after_KHI = (itKHI:Ns2D);
if 0
%% Phi frequency space time diagram at kz=0
fig = figure; FIGNAME = ['phi_freq_diag_',PARAMS];set(gcf, 'Position',  [100, 100, 500, 400]);
        [TY,TX] = meshgrid(Ts2D(after_KHI),kr(no_AA));
        pclr = pcolor(TX,TY,(squeeze(abs(PHI(no_AA,1,(after_KHI)))))); set(pclr, 'edgecolor','none'); colorbar;
        ylabel('$t c_s/R$'), xlabel('$0<k_r<2/3 k_r^{\max}$')
        legend('$\log|\tilde\phi(k_z=0)|$')
        title('Spectrogram of $\phi$')
end
%%
t0    = 0;
[~, it02D] = min(abs(Ts2D-t0));
[~, it05D] = min(abs(Ts5D-t0));
skip_ = 2; 
DELAY = 0.005*skip_;
FRAMES_2D = it02D:skip_:numel(Ts2D);
FRAMES_5D = it05D:skip_:numel(Ts5D);
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%% Density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 1;
FIELD = real(ni00); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$n_i$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Phi real space
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM),'_',PARAMS];INTERP = 1;
FIELD = real(phi); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$\phi$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% radial particle transport
GIFNAME = ['gamma_r',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 1;
FIELD = real(ni00.*dzphi); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$\Gamma_r$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Phi fourier
GIFNAME = ['FFT_phi',sprintf('_%.2d',JOBNUM),'_',PARAMS];INTERP = 0;
FIELD = ifftshift((abs(PHI)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$|\tilde\phi|$'; XNAME = '$k_r\rho_s$'; YNAME = '$k_z\rho_s$';
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
%% phi @ kz = 0
GIFNAME = ['phi_kz0',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0; SCALING = 0;
FIELD =squeeze(log10(abs(PHI(no_AA,1,:)))); linestyle = '-.'; FRAMES = FRAMES_2D;
X = kr(no_AA); T = Ts2D; YMIN = -30; YMAX = 6; XMIN = min(kr); XMAX = max(kr);
FIELDNAME = '$|\tilde\phi(k_z=0)|$'; XNAME = '$k_r\rho_s$';
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
FIELDNAME = '$N_i^{p0}$'; XNAME = '$k_{max}\rho_s$'; YNAME = '$P$';
create_gif_imagesc
end
if 0
%% maxkz, kr vs p, for all Nipj over time
GIFNAME = ['Nipj_kr',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
plt = @(x) squeeze(max((abs(x)),[],4));
FIELD = plt(Nipj); X = kr'; Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_i'; XNAME = '$k_r\rho_s$'; YNAME = '$P$, ${k_z}^{max}$';
create_gif_5D
end
if 0
%% maxkr, kz vs p, for all Nipj over time
GIFNAME = ['Nipj_kz',sprintf('_%.2d',JOBNUM),'_',PARAMS]; INTERP = 0;
plt = @(x) fftshift(squeeze(max((abs(x)),[],3)),3);
FIELD = plt(Nipj); X = sort(kz'); Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_i'; XNAME = '$k_z\rho_s$'; YNAME = '$P$, ${k_r}^{max}$';
create_gif_5D
end
%%