%% Load results
compile_results

%% Retrieving max polynomial degree and sampling info
Npe = numel(Pe); Nje = numel(Je); [JE,PE] = meshgrid(Je,Pe);
Npi = numel(Pi); Nji = numel(Ji); [JI,PI] = meshgrid(Ji,Pi);
Ns5D      = numel(Ts5D);
Ns2D      = numel(Ts2D);
% renaming and reshaping quantity of interest
Ts5D      = Ts5D';
Ts2D      = Ts2D';
if strcmp(OUTPUTS.write_non_lin,'.true.')
    Si00      = squeeze(Sipj(1,1,:,:,:));
    Se00      = squeeze(Sepj(1,1,:,:,:));
end
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

if strcmp(OUTPUTS.write_non_lin,'.true.')
for it = 1:numel(Ts5D)
    si00(:,:,it)      = real(fftshift(ifft2(squeeze(Si00(:,:,it)),Nr,Nz)));
end
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
if strcmp(OUTPUTS.write_non_lin,'.true.')
Se_norm  = zeros(Npe,Nje,Ns5D);% Time evol. of the norm of Sapj
Si_norm  = zeros(Npi,Nji,Ns5D);% .
Sne00_norm = zeros(1,Ns2D);    % Time evol. of the amp of e nonlin term
Sni00_norm = zeros(1,Ns2D);    %
end

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
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
if strcmp(OUTPUTS.write_non_lin,'.true.')   
    Se_norm(:,:,it)= sum(sum(abs(Sepj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Sne00_norm(it) = sum(sum(abs(Se00(:,:,it))))/Nkr/Nkz;
    Si_norm(:,:,it)= sum(sum(abs(Sipj(:,:,:,:,it)),3),4)/Nkr/Nkz;
    Sni00_norm(it) = sum(sum(abs(Si00(:,:,it))))/Nkr/Nkz;
end
end


%% Compute growth rate
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

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plots')
FMT = '.fig';
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM)];
    subplot(221); 
    for ip = 1:Npe
        for ij = 1:Nje
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_e^{',num2str(ip-1),num2str(ij-1),'}$'];
            semilogy(Ts5D,plt(Ne_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; xlabel('$t$'); ylabel('$\sum_{k_r,k_z}|N_e^{pj}|$');
    subplot(222)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            semilogy(Ts5D,plt(Ni_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; xlabel('$t$'); ylabel('$\sum_{k_r,k_z}|N_i^{pj}|$');
    subplot(223)
        plot(Ts2D,Flux_ri,'-','DisplayName','$\Gamma_{ri}$'); hold on;
        plot(Ts2D,Flux_zi,'-','DisplayName','$\Gamma_{zi}$'); hold on;
        plot(Ts2D,Flux_re,'-','DisplayName','$\Gamma_{re}$')
        plot(Ts2D,Flux_ze,'-','DisplayName','$\Gamma_{ze}$')
        grid on; xlabel('$t$'); ylabel('$\Gamma$'); %legend('show');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$S_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            semilogy(Ts5D,plt(Si_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; xlabel('$t$'); ylabel('$S$'); %legend('show');
else
%% Growth rate
    subplot(224)    
        [~,ikr0KH] = min(abs(kr-KR0KH));
        plot(kz(1:ikr0KH)/kr(ikr0KH),gkr0kz_Ni00(1:ikr0KH)/(KR0KH*A0KH),...
            'DisplayName',['J = ',num2str(JMAXI)]);
        xlabel('$k_z/k_{r0}$'); ylabel('$\gamma_{Ni}/(A_0k_{r0})$');
        xlim([0. 1.0]); %ylim([0. 0.04]);
end
save_figure

%%
if 1
%% Show frame in real space
tf = 100; [~,it] = min(abs(Ts2D-tf)); [~,it5D] = min(abs(Ts5D-tf));
fig = figure; FIGNAME = ['rz_frame',sprintf('_%.2d',JOBNUM)];
    subplot(221); plt = @(x) (((x)));
        pclr = pcolor((RR),(ZZ),plt(ne00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts2D(it))); legend('$|\hat n_e^{00}|$');
    subplot(222); plt = @(x) ((x));
        pclr = pcolor((RR),(ZZ),plt(ni00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts2D(it))); legend('$|\hat n_i^{00}|$');
    subplot(223); plt = @(x) ((x));
        pclr = pcolor((RR),(ZZ),plt(phi(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts2D(it))); legend('$|\hat\phi|$');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor((RR),(ZZ),plt(si00(:,:,it5D))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts5D(it5D))); legend('$|S_i^{00}|$');
end
save_figure
end

%%
t0    = 0;
skip_ = 1; 
DELAY = 0.01*skip_;
FRAMES = floor(t0/(Ts2D(2)-Ts2D(1)))+1:skip_:numel(Ts2D);
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%% Density ion
GIFNAME = ['ni00',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD = real(ni00); X = RR; Y = ZZ; T = Ts2D;
FIELDNAME = '$n_i^{00}$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
end
if 0
%% Density electron
GIFNAME = ['ne00',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ne00); X = RR; Y = ZZ; T = Ts2D;
FIELDNAME = '$n_e^{00}$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
end
if 0
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)];INTERP = 1;
FIELD = real(phi); X = RR; Y = ZZ; T = Ts2D;
FIELDNAME = '$\phi$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
end
if 0
%% Density ion frequency
GIFNAME = ['Ni00',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =ifftshift((abs(Ni00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts2D;
FIELDNAME = '$N_i^{00}$'; XNAME = '$k_r$'; YNAME = '$k_z$';
create_gif
end
if 0
%% Density ion frequency @ kr = 0
GIFNAME = ['Ni00_kr0',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =(squeeze(abs(Ni00(1,:,:)))); linestyle = 'o-.';
X = (kz); T = Ts2D; YMIN = -.1; YMAX = 1.1; XMIN = min(kz); XMAX = max(kz);
FIELDNAME = '$N_i^{00}(kr=0)$'; XNAME = '$k_r$';
create_gif_1D
% %% PJ ion moment frequency
% p_ = 0+1; j_ = 3+1;
% GIFNAME = ['Ni',num2str(p_),num2str(j_),sprintf('_%.2d',JOBNUM)]; INTERP = 0;
% FIELD =ifftshift(abs(squeeze(Nipj(p_,j_,:,:,:))),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts;
% FIELDNAME = ['$N_i^{',num2str(p_-1),num2str(j_-1),'}$']; XNAME = '$k_r$'; YNAME = '$k_z$';
% create_gif
% %% non linear term frequ. space
% GIFNAME = ['Si00',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
% FIELD = fftshift((abs(Si00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts;
% FIELDNAME = '$S_i^{00}$'; XNAME = '$k_r$'; YNAME = '$k_z$';
% create_gif
% %% Electron mode growth
% GIFNAME = ['norm_Nepj',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
% FIELD = Ne_norm; X = PE; Y = JE; T = Ts;
% FIELDNAME = '$\max_k{\gamma_e}$'; XNAME = '$p$'; YNAME = '$j$';
% create_gif
end
%%

if 0
%% Asymmetry error on kr = 0
up   = 2:Nkz/2;
down = Nkz:-1:Nkz/2+2;
eps_asym_r = Ts5D;
eps_asym_i = Ts5D;
for it = 1:numel(Ts5D)
    eps_asym_r(it) = max(squeeze(abs(real(Ni00(1,up,it))-real(Ni00(1,down,it)))));
    eps_asym_i(it) = max(abs(imag(Ni00(1,up,it))+imag(Ni00(1,down,it))));
end

figure
    it = 2;
    subplot(311)
        plot(kz(up),real(Ni00(1,up,it))); hold on
        plot(kz(up),real(Ni00(1,down,it)),'--')
        ylabel('Real'); title(['$N_i^{00}(t=',num2str(Ts5D(it)),')$'])
    subplot(312)
        plot(kz(up),imag(Ni00(1,up,it))); hold on
        plot(kz(up),-imag(Ni00(1,down,it)),'--')
        ylabel('Imag'); xlabel('$kz,kr=0$')
    subplot(313)
        plot(kz(2:Nkz/2),real(Ni00(1,up,it))-real(Ni00(1,down,it))); hold on
        plot(kz(2:Nkz/2),imag(Ni00(1,up,it))+imag(Ni00(1,down,it)),'--')
        legend('Re','Im');
        ylabel('err'); xlabel('$kz,kr=0$')
%%       
figure
     semilogy(Ts5D,eps_asym_r); hold on
     semilogy(Ts5D,eps_asym_i,'--');
     title('$\max_{k_z}|N_i^{00}(0,k_z)-N_i^{00}(0,-k_z)|$');
     legend('Re','Im'); grid on; xlabel('$t$')

end


if 0
%% Growth rate
fig = figure; FIGNAME = ['growth_rate',sprintf('_%.2d',JOBNUM)];
    [~,ikr0KH] = min(abs(kr-KR0KH));
    plot(kz/kr(ikr0KH),gkr0kz_Ni00/(kr(ikr0KH)*A0KH),'DisplayName',['J = ',num2str(JMAXI)])
    xlabel('$k_z/k_{r0}$'); ylabel('$\gamma_{Ni}/(A_0k_{r0})$');
    xlim([0. 1.0]); ylim([0. 0.6]);
save_figure

end
%%
if 0
%% Mode time evolution
[~,ik00] = min(abs(kz));
[~,idk]  = min(abs(kz-dkz));
[~,ik50] = min(abs(kz-0.5*KR0KH));
[~,ik75] = min(abs(kz-0.7*KR0KH));
[~,ik10] = min(abs(kz-1.00*KR0KH));
plt = @(x) abs(squeeze(x));
fig = figure; FIGNAME = ['mode_time_evolution',sprintf('_%.2d',JOBNUM)];
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik00,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik00)/KR0KH)]); hold on
        semilogy(Ts2D,plt(Ni00(ikr0KH,idk,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(idk)/KR0KH)]); hold on
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik50,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik50)/KR0KH)]); hold on
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik75,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik75)/KR0KH)]); hold on
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik10,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik10)/KR0KH)]); hold on
        xlabel('$t$'); ylabel('$\hat n_i^{00}$'); legend('show');
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik00,end))*exp(gkr0kz_Ni00(ik00)*(Ts2D-Ts2D(end))),'k--')
        semilogy(Ts2D,plt(Ni00(ikr0KH,idk,end))*exp(gkr0kz_Ni00(idk)*(Ts2D-Ts2D(end))),'k--')
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik50,end))*exp(gkr0kz_Ni00(ik50)*(Ts2D-Ts2D(end))),'k--')
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik75,end))*exp(gkr0kz_Ni00(ik75)*(Ts2D-Ts2D(end))),'k--')
        semilogy(Ts2D,plt(Ni00(ikr0KH,ik10,end))*exp(gkr0kz_Ni00(ik10)*(Ts2D-Ts2D(end))),'k--')
        plot(tstart*[1 1],ylim,'k-','LineWidth',0.5);
        plot(tend*[1 1],ylim,'k-','LineWidth',0.5);
save_figure
end

%%
if 0
%% Show frame in kspace
tf = 20; [~,it] = min(abs(Ts5D-tf));
fig = figure; FIGNAME = ['krkz_frame',sprintf('_%.2d',JOBNUM)];
    subplot(221); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(PHI(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts5D(it))); legend('$|\hat\phi|$');
    subplot(222); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ni00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts5D(it))); legend('$|\hat n_i^{00}|$');
    subplot(223); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ne00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts5D(it))); legend('$|\hat n_e^{00}|$');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Si00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_i^{00}$');
end
save_figure
end
