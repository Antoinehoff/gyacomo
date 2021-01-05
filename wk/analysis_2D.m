%% Load results
if 1
    %%
    outfile = '/marconi_scratch/userexternal/ahoffman/HeLaZ/results/Marconi/200x100_L_100_Pe_2_Je_1_Pi_2_Ji_1_nB_0.66_nN_1_nu_1e-01_FC_mu_1e-03/out.txt';
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
ne00   = zeros(Nr,Nz,Ns2D);
ni00   = zeros(Nr,Nz,Ns2D);
si00   = zeros(Nr,Nz,Ns5D);
phi    = zeros(Nr,Nz,Ns2D);
drphi  = zeros(Nr,Nz,Ns2D);
dzphi  = zeros(Nr,Nz,Ns2D);

for it = 1:numel(Ts2D)
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    ne00(:,:,it)  = real(fftshift(ifft2((NE_),Nr,Nz)));
    ni00(:,:,it)  = real(fftshift(ifft2((NI_),Nr,Nz)));
    phi (:,:,it)  = real(fftshift(ifft2((PH_),Nr,Nz)));
    drphi(:,:,it) = real(fftshift(ifft2(1i*KR.*(PH_),Nr,Nz)));
    dzphi(:,:,it) = real(fftshift(ifft2(1i*KZ.*(PH_),Nr,Nz)));
end

for it = 1:numel(Ts5D)
    si00(:,:,it)      = real(fftshift(ifft2(squeeze(Si00(:,:,it)),Nr,Nz)));
end

% Post processing
disp('- post processing')
E_pot    = zeros(1,Ns2D);      % Potential energy n^2
E_kin    = zeros(1,Ns2D);      % Kinetic energy grad(phi)^2
ExB      = zeros(1,Ns2D);      % ExB drift intensity \propto |\grad \phi|
Flux_ri  = zeros(1,Ns2D);      % Particle flux Gamma = <ni drphi>
Flux_zi  = zeros(1,Ns2D);      % Particle flux Gamma = <ni dzphi>
Flux_re  = zeros(1,Ns2D);      % Particle flux Gamma = <ne drphi>
Flux_ze  = zeros(1,Ns2D);      % Particle flux Gamma = <ne dzphi>
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
    Flux_ri(it)  = sum(sum(ni00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    Flux_zi(it)  = sum(sum(-ni00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
    Flux_re(it)  = sum(sum(ne00(:,:,it).*dzphi(:,:,it)))*dr*dz/Lr/Lz;
    Flux_ze(it)  = sum(sum(-ne00(:,:,it).*drphi(:,:,it)))*dr*dz/Lr/Lz;
end

E_kin_KZ = mean(mean(abs(Ddr.*PHI(:,:,it)).^2+abs(Ddz.*PHI(:,:,it)).^2,3),2);
E_kin_KR = mean(mean(abs(Ddr.*PHI(:,:,it)).^2+abs(Ddz.*PHI(:,:,it)).^2,3),2);
dEdt     = diff(E_pot+E_kin)./dt2D;

for it = 1:numel(Ts5D) % Loop over 5D arrays
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Se_norm(:,:,it)= sum(sum(abs(Sepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Si_norm(:,:,it)= sum(sum(abs(Sipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Sne00_norm(it) = sum(sum(abs(Se00(:,:,it))))/Nkr/Nkz;
    Sni00_norm(it) = sum(sum(abs(Si00(:,:,it))))/Nkr/Nkz;
end

%% Compute growth rate
if NON_LIN == 0
disp('- growth rate')
tend   = Ts2D(end); tstart   = 0.6*tend; 
g_          = zeros(Nkr,Nkz);
[~,ikr0KH] = min(abs(kr-KR0KH));
for ikr = 1:Nkr
    for ikz = 1:Nkz
        g_(ikr,ikz) = LinearFit_s(Ts2D,squeeze(abs(Ni00(ikr,ikz,:))),tstart,tend);
    end
end
% gkr0kz_Ni00 = max(real(g_(:,:)),[],1);
gkr0kz_Ni00 = real(g_(ikr0KH,:));
end
%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 0
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM)];
set(gcf, 'Position',  [100, 100, 900, 800])
    subplot(221); 
    for ip = 1:Npe
        for ij = 1:Nje
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_e^{',num2str(Pe(ip)),num2str(Je(ij)),'}$'];
            semilogy(Ts5D,plt(Ne_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_e^{pj}|$');
    subplot(222)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_i^{',num2str(Pi(ip)),num2str(Ji(ij)),'}$'];
            semilogy(Ts5D,plt(Ni_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_i^{pj}|$');
    subplot(223)
        plot(Ts2D,Flux_ri,'-','DisplayName','$\Gamma_{ri}$'); hold on;
        plot(Ts2D,Flux_zi,'-','DisplayName','$\Gamma_{zi}$'); hold on;
        plot(Ts2D,Flux_re,'-','DisplayName','$\Gamma_{re}$')
        plot(Ts2D,Flux_ze,'-','DisplayName','$\Gamma_{ze}$')
        grid on; xlabel('$t c_s/R$'); ylabel('$\Gamma$'); %legend('show');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$S_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            semilogy(Ts5D,plt(Si_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; xlabel('$t c_s/R$'); ylabel('$\sum_{k_r,k_z}|S_i^{pj}|$'); %legend('show');
else
%% Growth rate
    subplot(224)    
        [~,ikr0KH] = min(abs(kr-KR0KH));
        plot(kz(1:Nz/2),gkr0kz_Ni00(1:Nz/2),...
            'DisplayName',['J = ',num2str(JMAXI)]);
        xlabel('$k_z$'); ylabel('$\gamma_{Ni}$');
        xlim([0. 1.0]); %ylim([0. 0.04]);
end
suptitle(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB)]);
save_figure
end

%%
if 0
%% Photomaton : real space
tf = 0; [~,it] = min(abs(Ts2D-tf)); [~,it5D] = min(abs(Ts5D-tf));
fig = figure; FIGNAME = ['photo_real',sprintf('_t=%.0f',Ts2D(it))]; set(gcf, 'Position',  [100, 100, 1500, 500])
    subplot(131); plt = @(x) (((x))); 
        pclr = pcolor((RR),(ZZ),plt(ni00(:,:,it))); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        xlabel('$r/\rho_s$'); ylabel('$z/\rho_s$'); legend('$n_i$');
        
    subplot(132); plt = @(x) ((x));
        DATA = plt(ni00(:,:,it))-plt(ne00(:,:,it));
        pclr = pcolor((RR),(ZZ),DATA./max(max(DATA))); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        xlabel('$r/\rho_s$'); legend('$n_i-n_e$'); set(gca,'ytick',[]);
        
        
    subplot(133); plt = @(x) ((x));
        DATA = plt(phi(:,:,it));
        pclr = pcolor((RR),(ZZ),DATA./max(max(DATA))); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        xlabel('$r/\rho_s$'); set(gca,'ytick',[]); legend('$\phi$');
        
% if strcmp(OUTPUTS.write_non_lin,'.true.')
%     subplot(133); plt = @(x) fftshift((abs(x)),2);
%         pclr = pcolor((RR),(ZZ),plt(si00(:,:,it5D))); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
%         xlabel('$r/\rho_s$'); legend('$|S_i^{00}|$'); set(gca,'ytick',[])
% end
suptitle(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB), sprintf(', $t c_s/R=%.0f$',Ts2D(it))]);
save_figure
end

if 0
%% Photomaton : real space
FIELD = ni00; FNAME = 'ni';
% FIELD = ne00; FNAME = 'ne';
% FIELD = phi; FNAME = 'phi';
tf = 19;  [~,it1] = min(abs(Ts2D-tf));
tf = 20;  [~,it2] = min(abs(Ts2D-tf)); 
tf = 21; [~,it3] = min(abs(Ts2D-tf));
tf = 22; [~,it4] = min(abs(Ts2D-tf));
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
t0    = 40;
skip_ = 1; 
DELAY = 0.01*skip_;
FRAMES = floor(t0/(Ts2D(2)-Ts2D(1)))+1:skip_:numel(Ts2D);
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%% Density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ni00); X = RR; Y = ZZ; T = Ts2D;
FIELDNAME = '$n_i$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Density electron
GIFNAME = ['ne',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ne00); X = RR; Y = ZZ; T = Ts2D;
FIELDNAME = '$n_e$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Density ion - electron
GIFNAME = ['ni-ne',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ni00+ne00); X = RR; Y = ZZ; T = Ts2D;
FIELDNAME = '$n_i-n_e$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)];INTERP = 1;
FIELD = real(phi); X = RR; Y = ZZ; T = Ts2D;
FIELDNAME = '$\phi$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end
if 0
%% Density ion frequency
GIFNAME = ['Ni00',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =ifftshift((abs(Ni00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts2D;
FIELDNAME = '$N_i^{00}$'; XNAME = '$k_r\rho_s$'; YNAME = '$k_z\rho_s$';
create_gif
end
if 0
%% Density ion frequency @ kr = 0
GIFNAME = ['Ni00_kr0',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =(squeeze(abs(Ni00(1,:,:)))); linestyle = 'o-.';
X = (kz); T = Ts2D; YMIN = -.1; YMAX = 1.1; XMIN = min(kz); XMAX = max(kz);
FIELDNAME = '$N_i^{00}(kr=0)$'; XNAME = '$k_r\rho_s$';
create_gif_1D
end
%%

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
t0  = 1; t1 = 400; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
fig = figure; FIGNAME = 'space_time_drphi';set(gcf, 'Position',  [100, 100, 1200, 400])
    subplot(211)
    plot(Ts2D(it0:it1),Flux_ri(it0:it1));
    ylabel('$\Gamma_r$'); grid on
%     title(['$\eta=',num2str(ETAB),'\quad',...
%         '\nu_{',CONAME,'}=',num2str(NU),'$'])
%     legend(['$P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$'])
    ylim([0,1.1*max(Flux_ri(it0:it1))]);
    subplot(212)
    [TY,TX] = meshgrid(r,Ts2D(it0:it1));
    pclr = pcolor(TX,TY,squeeze(mean(drphi(:,:,it0:it1),2))'); set(pclr, 'edgecolor','none'); %colorbar;
    xlabel('$t c_s/R$'), ylabel('$r/\rho_s$')
    legend('$\langle\partial_r \phi\rangle_z$')
save_figure
end

%%
if 0
%% Mode time evolution
[~,ikr ] = min(abs(kr-dkr));
[~,ik00] = min(abs(kz));
[~,idk]  = min(abs(kz-dkz));
[~,ik50] = min(abs(kz-0.1*max(kz)));
[~,ik75] = min(abs(kz-0.2*max(kz)));
[~,ik10] = min(abs(kz-0.3*max(kz)));
plt = @(x) abs(squeeze(x));
fig = figure; FIGNAME = ['mode_time_evolution',sprintf('_%.2d',JOBNUM)];
        semilogy(Ts2D,plt(Ni00(ikr,ik00,:)),'DisplayName', ...
            ['$k_z = $',num2str(kz(ik00))]); hold on
        semilogy(Ts2D,plt(Ni00(ikr,idk,:)),'DisplayName', ...
            ['$k_z = $',num2str(kz(idk))]); hold on
        semilogy(Ts2D,plt(Ni00(ikr,ik50,:)),'DisplayName', ...
            ['$k_z = $',num2str(kz(ik50))]); hold on
        semilogy(Ts2D,plt(Ni00(ikr,ik75,:)),'DisplayName', ...
            ['$k_z = $',num2str(kz(ik75))]); hold on
        semilogy(Ts2D,plt(Ni00(ikr,ik10,:)),'DisplayName', ...
            ['$k_z = $',num2str(kz(ik10))]); hold on
        xlabel('$t$'); ylabel('$\hat n_i^{00}$'); legend('show');
title(sprintf('$k_r=$ %1.1f',kr(ikr)))
save_figure
end

%%
if 0
%% Show frame in kspace
tf = 300; [~,it2] = min(abs(Ts2D-tf)); [~,it5] = min(abs(Ts5D-tf));
fig = figure; FIGNAME = ['krkz_frame',sprintf('t=%.0f',Ts2D(it2))];set(gcf, 'Position',  [100, 100, 700, 600])
    subplot(221); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(PHI(:,:,it2))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('$t c_s/R=%.0f$',Ts2D(it2))); legend('$|\hat\phi|$');
    subplot(222); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ni00(:,:,it2))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_i^{00}|$');
    subplot(223); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ne00(:,:,it2))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_e^{00}|$');
    subplot(224); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Si00(:,:,it5))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_i^{00}$');
save_figure
end

%% Phase space distribution function
% M_ = 25;
% spar = linspace(0,4,M_); xperp = spar;
% 
% PSDF_e = zeros(numel(kr),numel(kz),M_,M_,numel(Ts5D));
% for ikr = 1:numel(kr)
%     for ikz = 1:numel(kz)
%         for it = 1:numel(Ts5D)
%         PSDF_e(ikr,ikz,:,:,it) = compute_fa(Nepj(:,:,ikr,ikz,it), spar, xperp);
%         end
%     end
% end
% 
% ktarget = 0.3*max(kz);
% ttarget = 310;
% 
% [~,ik10] = min(abs(kz-ktarget));
% [~,it]   = min(abs(Ts5D-ttarget));
% if 0
%     %%
% plt = @(x) real(x(ik10,ik10,:,:,it));
% pclr = pcolor(spar,xperp,plt(PSDF_e)); set(pclr, 'edgecolor','none'); colorbar;
% xlabel('$s_\parallel$'); ylabel('$x_\perp$'); title(sprintf('$t c_s/R=%.0f$',Ts5D(it))); 
% legend(['$Re(f_e),k_\perp \approx$',sprintf('%01.0f',norm([kr(ik10),kz(ik10)]))]);
% end
% if 0
% %% Phase space distribution function time evolution
% GIFNAME = ['f_e',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
% FIELD = real(ni00+ne00); X = RR; Y = ZZ; T = Ts5D;
% FIELDNAME = '$n_i-n_e$'; XNAME = '$r\rho_s$'; YNAME = '$z\rho_s$';
% create_gif
% end

