%% Load results
if 0
    %%
    outfile ='';
    outfile ='';
    outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/Marconi_DGGK_eta_0.6_nu_1e+00/150x75_L_70_P_10_J_5_eta_0.6_nu_1e+00_DGGK_CLOS_0_mu_8e-04/out.txt';

    BASIC.RESDIR = load_marconi(outfile);
end
%%
% JOBNUM = 0; load_results;
% JOBNUM = 1; load_results;
compile_results
load_params

%% Retrieving max polynomial degree and sampling info
Npe = numel(Pe); Nje = numel(Je); [JE,PE] = meshgrid(Je,Pe);
Npi = numel(Pi); Nji = numel(Ji); [JI,PI] = meshgrid(Ji,Pi);
Ns5D      = numel(Ts5D);
Ns2D      = numel(Ts2D);
% renaming and reshaping quantity of interest
Ts5D      = Ts5D';
Ts2D      = Ts2D';
Si00      = squeeze(Sipj(1,1,:,:,:));
Se00      = squeeze(Sepj(1,1,:,:,:));
%% Build grids
Nkr = numel(kr); Nkz = numel(kz);
[KZ,KR] = meshgrid(kz,kr);
Lkr = max(kr)-min(kr); Lkz = max(kz)-min(kz);
dkr = Lkr/(Nkr-1); dkz = Lkz/(Nkz-1);
KPERP2 = KZ.^2+KR.^2;

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

for it = 1:numel(Ts5D)
    [~, it2D] = min(abs(Ts2D-Ts5D(it)));
    si00(:,:,it)      = real(fftshift(ifft2(squeeze(Si00(:,:,it)),Nr,Nz)));
    
    Np_i = zeros(Nkr,Nkz); % Ion particle density in Fourier space
    for ij = 1:Nji
        Kn = (KPERP2/2.).^(ij-1) .* exp(-KPERP2/2)/(factorial(ij-1));
        Np_i = Np_i + Kn.*squeeze(Nipj(1,ij,:,:,it));
    end
    np_i(:,:,it)      = real(fftshift(ifft2(squeeze(Np_i(:,:)),Nr,Nz)));
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

Ne_norm  = zeros(Npe,Nje,Ns5D);% Time evol. of the norm of Napj
Ni_norm  = zeros(Npi,Nji,Ns5D);% .
Se_norm  = zeros(Npe,Nje,Ns5D);% Time evol. of the norm of Sapj
Si_norm  = zeros(Npi,Nji,Ns5D);% .
Sne00_norm = zeros(1,Ns2D);    % Time evol. of the amp of e nonlin term
Sni00_norm = zeros(1,Ns2D);    %


Ddr = 1i*KR; Ddz = 1i*KZ; lapl   = Ddr.^2 + Ddz.^2; 

for it = 1:numel(Ts2D) % Loop over 2D arrays
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    E_pot(it)   = pi/Lr/Lz*sum(sum(abs(NI_).^2))/Nkr/Nkr; % integrate through Parseval id
    E_kin(it)   = pi/Lr/Lz*sum(sum(abs(Ddr.*PH_).^2+abs(Ddz.*PH_).^2))/Nkr/Nkr;
    ExB(it)     = max(max(max(abs(phi(3:end,:,it)-phi(1:end-2,:,it))/(2*dr))),max(max(abs(phi(:,3:end,it)-phi(:,1:end-2,it))'/(2*dz))));
    GFlux_ri(it)  = sum(sum(ni00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_zi(it)  = sum(sum(-ni00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_re(it)  = sum(sum(ne00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    GFlux_ze(it)  = sum(sum(-ne00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
end

E_kin_KZ = mean(mean(abs(Ddr.*PHI(:,:,it)).^2+abs(Ddz.*PHI(:,:,it)).^2,3),2);
E_kin_KR = mean(mean(abs(Ddr.*PHI(:,:,it)).^2+abs(Ddz.*PHI(:,:,it)).^2,3),2);
dEdt     = diff(E_pot+E_kin)./dt2D;

for it = 1:numel(Ts5D) % Loop over 5D arrays
    [~, it2D] = min(abs(Ts2D-Ts5D(it)));
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Se_norm(:,:,it)= sum(sum(abs(Sepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Si_norm(:,:,it)= sum(sum(abs(Sipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Sne00_norm(it) = sum(sum(abs(Se00(:,:,it))))/Nkr/Nkz;
    Sni00_norm(it) = sum(sum(abs(Si00(:,:,it))))/Nkr/Nkz;
    % Particle flux
    PFlux_ri(it)   = sum(sum(np_i(:,:,it).*dzphi(:,:,it2D)))*dr*dz/Lr/Lz;
end

%% Compute growth rate
disp('- growth rate')
% Find max value of transport (end of linear mode)
[~,itmax] = max(GFlux_ri);
tstart = 0.1 * Ts2D(itmax); tend = 0.9 * Ts2D(itmax);
g_          = zeros(Nkr,Nkz);
[~,ikr0KH] = min(abs(kr-KR0KH));
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
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM)];
set(gcf, 'Position',  [100, 100, 900, 800])
    subplot(221); 
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
    subplot(222)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_i^{',num2str(Pi(ip)),num2str(Ji(ij)),'}$'];
            clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
            lstyle   = line_styles(min(ij,numel(line_styles)));
            plot(Ts5D,plt(Ni_norm),'DisplayName',plotname,...
                'Color',clr,'LineStyle',lstyle{1}); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_i^{pj}|$');
    subplot(223)
        plot(kz,g_(1,:),'-','DisplayName','$\gamma$'); hold on;
        grid on; xlabel('$k_z\rho_s$'); ylabel('$\gamma R/c_s$'); %legend('show');
    subplot(224)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$S_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
            lstyle   = line_styles(min(ij,numel(line_styles)));
            semilogy(Ts5D,plt(Si_norm),'DisplayName',plotname,...
                'Color',clr,'LineStyle',lstyle{1}); hold on;
        end
    end
    grid on; xlabel('$t c_s/R$'); ylabel('$\sum_{k_r,k_z}|S_i^{pj}|$'); %legend('show');
% suptitle(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB)]);
save_figure
end

if 1
%% Particle fluxes
fig = figure; FIGNAME = ['gamma',sprintf('_%.2d',JOBNUM)];
set(gcf, 'Position',  [100, 100, 1200, 400])
    subplot(211)
        plot(Ts2D,GFlux_ri); hold on
        plot(Ts5D,PFlux_ri,'--'); hold on
        ylabel('$\Gamma_r$'); grid on
        title(['$\eta=',num2str(ETAB),'\quad',...
            '\nu_{',CONAME,'}=',num2str(NU),'$'])
        legend(['$P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$'],'Particle flux')%'$\eta\gamma_{max}/k_{max}^2$')
        set(gca,'xticklabel',[])
    subplot(212)
        plot(Ts2D,GFLUX_RI); hold on
        plot(Ts5D,PFLUX_RI,'--'); hold on
        ylabel('$\Gamma_r$'); grid on
        title(['$\eta=',num2str(ETAB),'\quad',...
            '\nu_{',CONAME,'}=',num2str(NU),'$'])
        legend(['$P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$'],'Particle flux')%'$\eta\gamma_{max}/k_{max}^2$')
        set(gca,'xticklabel',[])
save_figure
end

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
fig = figure; FIGNAME = 'space_time_drphi';set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
        plot(Ts2D,GFlux_ri); hold on
        plot(Ts5D,PFlux_ri,'.'); hold on
        plot(Ts2D,GFLUX_RI,'--'); hold on
%         plot(Ts2D,Bohm_transport*ones(size(Ts2D)),'--'); hold on
        ylabel('$\Gamma_r$'); grid on
        title(['$\eta=',num2str(ETAB),'\quad',...
            '\nu_{',CONAME,'}=',num2str(NU),'$'])
        legend(['$P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$'],'Particle flux')%'$\eta\gamma_{max}/k_{max}^2$')
        set(gca,'xticklabel',[])
    subplot(312)
        yyaxis left
        plot(Ts2D,squeeze(max(max((phi)))))
        ylabel('$\max \phi$')
        yyaxis right
        plot(Ts2D,squeeze(mean(max(dr2phi))))
        ylabel('$s\sim\langle\partial_r^2\phi\rangle_z$'); grid on  
        set(gca,'xticklabel',[])
    subplot(313)
        [TY,TX] = meshgrid(r,Ts2D);
        pclr = pcolor(TX,TY,squeeze(mean(drphi(:,:,:),2))'); set(pclr, 'edgecolor','none'); %colorbar;
        xlabel('$t c_s/R$'), ylabel('$r/\rho_s$')
        legend('$\langle\partial_r \phi\rangle_z$')
save_figure
end

if 0
%% Photomaton : real space
% FIELD = ni00; FNAME = 'ni';
% FIELD = ne00; FNAME = 'ne';
FIELD = phi; FNAME = 'phi';
tf = 200;  [~,it1] = min(abs(Ts2D-tf));
tf = 600;  [~,it2] = min(abs(Ts2D-tf)); 
tf =1000; [~,it3] = min(abs(Ts2D-tf));
tf =2000; [~,it4] = min(abs(Ts2D-tf));
fig = figure; FIGNAME = [FNAME,'_snaps']; set(gcf, 'Position',  [100, 100, 1500, 400])
plt = @(x) x;%./max(max(x));
    subplot(141)
        DATA = plt(FIELD(:,:,it1));
        pclr = pcolor((RR),(ZZ),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        xlabel('$r/\rho_s$'); ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it1)));
    subplot(142)
        DATA = plt(FIELD(:,:,it2));
        pclr = pcolor((RR),(ZZ),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$'); set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it2)));
    subplot(143)
        DATA = plt(FIELD(:,:,it3));
        pclr = pcolor((RR),(ZZ),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it3)));
    subplot(144)
        DATA = plt(FIELD(:,:,it4));
        pclr = pcolor((RR),(ZZ),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
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
fig = figure; FIGNAME = ['krkz_',sprintf('t=%.0f',Ts2D(it2))];set(gcf, 'Position',  [100, 100, 700, 600])
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
if 1
%% Ion moments max mode vs pj
% tf = Ts2D(end-3); 
for tf = []
[~,it2] = min(abs(Ts2D-tf)); [~,it5] = min(abs(Ts5D-tf));
% it2 = it2 + 1;
fig = figure; FIGNAME = ['kmaxp_Nipj_',sprintf('t=%.2f',Ts2D(it2))];set(gcf, 'Position',  [100, 100, 700, 600])

plt = @(x) squeeze(max(abs(x),[],4));
% plt = @(x) squeeze(max(fftshift(abs(x),2),[],4));

for ij_ = 1:numel(Ji)
    subplot(100+numel(Ji)*10+ij_)
        pclr = imagesc(kr,Pi,plt(Nipj(:,ij_,:,:,it5)));
        xlabel('$k_r$');
        if ij_ == 1
            ylabel('$P$(max o. $k_z$)');
        else
            yticks([])
        end
        LEGEND = ['$|\hat n_i^{p',num2str(ij_-1),'}|$']; title(LEGEND);
end
save_figure
end
end


%%
t0    = 0;
[~, it02D] = min(abs(Ts2D-t0));
[~, it05D] = min(abs(Ts5D-t0));
skip_ = 1; 
DELAY = 0.02*skip_;
FRAMES_2D = it02D:skip_:numel(Ts2D);
FRAMES_5D = it05D:skip_:numel(Ts5D);
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%% Density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ni00); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$n_i$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Density electron
GIFNAME = ['ne',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ne00); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$n_e$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Phi real space
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)];INTERP = 1;
FIELD = real(phi); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$\phi$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Phi fourier
GIFNAME = ['FFT_phi',sprintf('_%.2d',JOBNUM)];INTERP = 0;
FIELD = ifftshift((abs(PHI)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$|\tilde\phi|$'; XNAME = '$k_r\rho_s$'; YNAME = '$k_z\rho_s$';
create_gif
end
if 0
%% phi @ z = 0
GIFNAME = ['phi_r0',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =(squeeze(real(phi(:,1,:)))); linestyle = '-.'; FRAMES = FRAMES_2D;
X = (r); T = Ts2D; YMIN = -1.1; YMAX = 1.1; XMIN = min(r); XMAX = max(r);
FIELDNAME = '$\phi(r=0)$'; XNAME = '$r/\rho_s$';
create_gif_1D
end
if 0
%% Density ion frequency
GIFNAME = ['Ni00',sprintf('_%.2d',JOBNUM)]; INTERP = 0; FRAMES = FRAMES_2D;
FIELD =ifftshift((abs(Ni00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts2D;
FIELDNAME = '$N_i^{00}$'; XNAME = '$k_r\rho_s$'; YNAME = '$k_z\rho_s$';
create_gif
end
if 0
%% Density ion frequency @ kr = 0
GIFNAME = ['Ni00_kr0',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =(squeeze(abs(Ni00(1,:,:)))); linestyle = 'o-.'; FRAMES = FRAMES_2D;
X = (kz); T = Ts2D; YMIN = -.1; YMAX = 1.1; XMIN = min(kz); XMAX = max(kz);
FIELDNAME = '$N_i^{00}(kr=0)$'; XNAME = '$k_r\rho_s$';
create_gif_1D
end
if 0
%% kr vs P Si
GIFNAME = ['Sip0_kr',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
plt = @(x) squeeze(max((abs(x)),[],4));
FIELD =plt(Sipj(:,1,:,:,:)); X = kr'; Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = '$N_i^{p0}$'; XNAME = '$k_{max}\rho_s$'; YNAME = '$P$';
create_gif_imagesc
end
if 1
%% maxkz, kr vs p, for all Nipj over time
GIFNAME = ['Nipj_kr',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
plt = @(x) squeeze(max((abs(x)),[],4));
FIELD = plt(Nipj); X = kr'; Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_i'; XNAME = '$k_r\rho_s$'; YNAME = '$P$, ${k_z}^{max}$';
create_gif_5D
end
if 1
%% maxkr, kz vs p, for all Nipj over time
GIFNAME = ['Nipj_kz',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
plt = @(x) fftshift(squeeze(max((abs(x)),[],3)),3);
FIELD = plt(Nipj); X = sort(kz'); Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_i'; XNAME = '$k_z\rho_s$'; YNAME = '$P$, ${k_r}^{max}$';
create_gif_5D
end
if 0
%% maxkz, kr vs p, for all Nepj over time
GIFNAME = ['Nepj_kr',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
plt = @(x) squeeze(max((abs(x)),[],4));
FIELD = plt(Nepj); X = kr'; Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_e'; XNAME = '$k_r\rho_s$'; YNAME = '$P$, ${k_z}^{max}$';
create_gif_5D
end
if 0
%% maxkz, kz vs p, for all Nepj over time
GIFNAME = ['Nepj_kz',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
plt = @(x) fftshift(squeeze(max((abs(x)),[],3)),3);
FIELD = plt(Nepj); X = sort(kz'); Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_e'; XNAME = '$k_z\rho_s$'; YNAME = '$P$, ${k_r}^{max}$';
create_gif_5D
end
%%