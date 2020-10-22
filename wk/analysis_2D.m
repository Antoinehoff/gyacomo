%% load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JOBNUM = 00;
FMT = '.fig';
filename = [BASIC.SIMID,'_','%.2d.h5'];
filename = sprintf(filename,JOBNUM); disp(['Loading ',filename])
% Loading from output file
if strcmp(OUTPUTS.write_moments,'.true.') 
[Nipj, Pi, Ji, kr, kz, Ts, dt] = load_5D_data(filename, 'moments_i');
[Nepj, Pe, Je,  ~,  ~,  ~,  ~] = load_5D_data(filename, 'moments_e');
Ni00    = squeeze(Nipj(1,1,:,:,:));
Ne00    = squeeze(Nepj(1,1,:,:,:));
else
[Ni00, kr, kz, Ts, dt] = load_2D_data(filename, 'Ni00');
 Ne00                  = load_2D_data(filename, 'Ne00');
 Pi = [0]; Ji = Pi; Pe = Pi; Je = Pi;
 Nipj = zeros(1,1,numel(kr),numel(kz),numel(Ts));
 Nepj = Nipj;
 Nipj(1,1,:,:,:) = Ni00; Nepj(1,1,:,:,:) = Ne00;
end
PHI                            = load_2D_data(filename, 'phi');

if strcmp(OUTPUTS.write_non_lin,'.true.') 
    Sipj    = load_5D_data(filename, 'Sipj');
    Sepj    = load_5D_data(filename, 'Sepj');
end
%% Retrieving max polynomial degree and sampling info
Npe = numel(Pe); Nje = numel(Je); [JE,PE] = meshgrid(Je,Pe);
Npi = numel(Pi); Nji = numel(Ji); [JI,PI] = meshgrid(Ji,Pi);
dt_samp = mean(diff(Ts)); Ns      = numel(Ts);
% renaming and reshaping quantity of interest
Ts      = Ts';
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
ne00   = zeros(Nr,Nz);
ni00   = zeros(Nr,Nz);
si00   = zeros(Nr,Nz);
phi    = zeros(Nr,Nz);

for it = 1:numel(PHI(1,1,:))
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    ne00(:,:,it)  = real(fftshift(ifft2((NE_),Nr,Nz)));
    ni00(:,:,it)  = real(fftshift(ifft2((NI_),Nr,Nz)));
    phi (:,:,it)  = real(fftshift(ifft2((PH_),Nr,Nz)));
if strcmp(OUTPUTS.write_non_lin,'.true.')
    SI_ = Si00(:,:,it);
    si00(:,:,it)  = fftshift(ifft2(half_2_full_cc_2D(SI_),'symmetric'));
end
end

% Post processing
disp('- post processing')
E_pot    = zeros(1,Ns);      % Potential energy n^2
E_kin    = zeros(1,Ns);      % Kinetic energy grad(phi)^2
ExB      = zeros(1,Ns);      % ExB drift intensity \propto |\grad \phi|
CFL      = zeros(1,Ns);      % CFL
Ne_norm  = zeros(Npe,Nje,Ns);% Time evol. of the norm of Napj
Ni_norm  = zeros(Npi,Nji,Ns);% .
if strcmp(OUTPUTS.write_non_lin,'.true.')
Se_norm  = zeros(Npe,Nje,Ns);% Time evol. of the norm of Sapj
Si_norm  = zeros(Npi,Nji,Ns);% .
Sne00_norm = zeros(1,Ns);    % Time evol. of the amp of e nonlin term
Sni00_norm = zeros(1,Ns);    %
end

Ddr = 1i*KR; Ddz = 1i*KZ; lapl   = Ddr.^2 + Ddz.^2; 

for it = 1:numel(PHI(1,1,:))
    NE_ = Ne00(:,:,it); NI_ = Ni00(:,:,it); PH_ = PHI(:,:,it);
    Ne_norm(:,:,it)= sum(sum(abs(Nepj(:,:,:,:,it)),3),4);
    Ni_norm(:,:,it)= sum(sum(abs(Nipj(:,:,:,:,it)),3),4);
if strcmp(OUTPUTS.write_non_lin,'.true.')   
    Se_norm(:,:,it)= sum(sum(abs(Sepj(:,:,:,:,it)),3),4);
    Sne00_norm(it) = sum(sum(abs(Se00(:,:,it))));
    Si_norm(:,:,it)= sum(sum(abs(Sipj(:,:,:,:,it)),3),4);
    Sni00_norm(it) = sum(sum(abs(Si00(:,:,it))));
end
    E_pot(it)   = pi/Lr/Lz*sum(sum(abs(NI_).^2))/Nkr/Nkr; % integrate through Parseval id
    E_kin(it)   = pi/Lr/Lz*sum(sum(abs(Ddr.*PH_).^2+abs(Ddz.*PH_).^2))/Nkr/Nkr;
    ExB(it)     = max(max(max(abs(phi(3:end,:,it)-phi(1:end-2,:,it))/(2*dr))),max(max(abs(phi(:,3:end,it)-phi(:,1:end-2,it))'/(2*dz))));
end
E_kin_KZ = mean(mean(abs(Ddr.*PHI(:,:,it)).^2+abs(Ddz.*PHI(:,:,it)).^2,3),2);
E_kin_KR = mean(mean(abs(Ddr.*PHI(:,:,it)).^2+abs(Ddz.*PHI(:,:,it)).^2,3),2);
dEdt     = diff(E_pot+E_kin)./diff(Ts);

%% Growth rate
disp('- growth rate')
tend   = Ts(end); tstart   = 0.8*tend; 
g_          = zeros(Nkr,Nkz);
[~,ikr0KH] = min(abs(kr-KR0KH));
for ikr = 1:Nkr
    for ikz = 1:Nkz
        g_(ikr,ikz) = LinearFit_s(Ts,squeeze(abs(Ni00(ikr,ikz,:))),tstart,tend);
    end
end
% gkr0kz_Ni00 = max(real(g_(:,:)),[],1);
gkr0kz_Ni00 = real(g_(ikr0KH,:));

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plots')
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM)];
    subplot(221); 
    for ip = 1:Npe
        for ij = 1:Nje
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_e^{',num2str(ip-1),num2str(ij-1),'}$'];
            semilogy(Ts,plt(Ne_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; xlabel('$t$'); ylabel('$\sum_{k_r,k_z}|N_e^{pj}|$');
    subplot(222)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            semilogy(Ts,plt(Ni_norm),'DisplayName',plotname); hold on;
        end
    end
    grid on; xlabel('$t$'); ylabel('$\sum_{k_r,k_z}|N_i^{pj}|$');
    subplot(223)
        semilogy(Ts,E_kin+E_pot,'-','DisplayName','$\sum|ik\tilde\phi_i|^2+\sum|\tilde n_i|^2$')
        hold on;
        grid on; xlabel('$t$'); ylabel('$E$'); %legend('show');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224)
    for ip = 1:Npi
        for ij = 1:Nji
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$S_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            semilogy(Ts,plt(Si_norm),'DisplayName',plotname); hold on;
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
t0    = 0;
skip_ = 1; 
DELAY = 0.01*skip_;
FRAMES = floor(t0/dt_samp)+1:skip_:numel(Ts);
if 0
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Density ion
GIFNAME = ['ni00',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ni00); X = RR; Y = ZZ; T = Ts;
FIELDNAME = '$n_i^{00}$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
%% Density electron
GIFNAME = ['ne00',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(ne00); X = RR; Y = ZZ; T = Ts;
FIELDNAME = '$n_e^{00}$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)];INTERP = 1;
FIELD = real(phi); X = RR; Y = ZZ; T = Ts;
FIELDNAME = '$\phi$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
%% Density ion frequency
GIFNAME = ['Ni00',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =ifftshift((abs(Ni00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts;
FIELDNAME = '$N_i^{00}$'; XNAME = '$k_r$'; YNAME = '$k_z$';
create_gif
%% PJ ion moment frequency
p_ = 0+1; j_ = 1+1;
GIFNAME = ['Ni',num2str(p_),num2str(j_),sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD =ifftshift(abs(squeeze(Nipj(p_,j_,:,:,:))),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts;
FIELDNAME = ['$N_i^{',num2str(p_-1),num2str(j_-1),'}$']; XNAME = '$k_r$'; YNAME = '$k_z$';
create_gif
%% non linear term frequ. space
GIFNAME = ['Si00',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
FIELD = fftshift((abs(Si00)),2); X = fftshift(KR,2); Y = fftshift(KZ,2); T = Ts;
FIELDNAME = '$S_i^{00}$'; XNAME = '$k_r$'; YNAME = '$k_z$';
create_gif
%% non linear term real space
GIFNAME = ['si00',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = real(si00); X = RR; Y = ZZ; T = Ts;
FIELDNAME = '$s_i^{00}$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
%% Electron mode growth
GIFNAME = ['norm_Nepj',sprintf('_%.2d',JOBNUM)]; INTERP = 1;
FIELD = Ne_norm; X = PE; Y = JE; T = Ts;
FIELDNAME = '$\max_k{\gamma_e}$'; XNAME = '$p$'; YNAME = '$j$';
create_gif
end
%%

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
if 1
%% Mode time evolution
[~,ik00] = min(abs(kz));
[~,idk]  = min(abs(kz-dkz));
[~,ik50] = min(abs(kz-0.5*KR0KH));
[~,ik75] = min(abs(kz-0.7*KR0KH));
[~,ik10] = min(abs(kz-1.00*KR0KH));
plt = @(x) abs(squeeze(x));
fig = figure; FIGNAME = ['frame',sprintf('_%.2d',JOBNUM)];
        semilogy(Ts,plt(Ni00(ikr0KH,ik00,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik00)/KR0KH)]); hold on
        semilogy(Ts,plt(Ni00(ikr0KH,idk,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(idk)/KR0KH)]); hold on
        semilogy(Ts,plt(Ni00(ikr0KH,ik50,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik50)/KR0KH)]); hold on
        semilogy(Ts,plt(Ni00(ikr0KH,ik75,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik75)/KR0KH)]); hold on
        semilogy(Ts,plt(Ni00(ikr0KH,ik10,:)),'DisplayName', ...
            ['$k_z/k_{r0} = $',num2str(kz(ik10)/KR0KH)]); hold on
        xlabel('$t$'); ylabel('$\hat n_i^{00}$'); legend('show');
        semilogy(Ts,plt(Ni00(ikr0KH,ik00,end))*exp(gkr0kz_Ni00(ik00)*(Ts-Ts(end))),'k--')
        semilogy(Ts,plt(Ni00(ikr0KH,idk,end))*exp(gkr0kz_Ni00(idk)*(Ts-Ts(end))),'k--')
        semilogy(Ts,plt(Ni00(ikr0KH,ik50,end))*exp(gkr0kz_Ni00(ik50)*(Ts-Ts(end))),'k--')
        semilogy(Ts,plt(Ni00(ikr0KH,ik75,end))*exp(gkr0kz_Ni00(ik75)*(Ts-Ts(end))),'k--')
        semilogy(Ts,plt(Ni00(ikr0KH,ik10,end))*exp(gkr0kz_Ni00(ik10)*(Ts-Ts(end))),'k--')
        plot(tstart*[1 1],ylim,'k-','LineWidth',0.5);
        plot(tend*[1 1],ylim,'k-','LineWidth',0.5);
save_figure
end

%%
if 0
%% Show frame in real space
tf = 20; [~,it] = min(abs(Ts-tf));
fig = figure; FIGNAME = ['rz_frame',sprintf('_%.2d',JOBNUM)];
    subplot(221); plt = @(x) (((x)));
        pclr = pcolor((RR),(ZZ),plt(ne00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts(it))); legend('$|\hat n_e^{00}|$');
    subplot(222); plt = @(x) ((x));
        pclr = pcolor((RR),(ZZ),plt(ni00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts(it))); legend('$|\hat n_i^{00}|$');
    subplot(223); plt = @(x) ((x));
        pclr = pcolor((RR),(ZZ),plt(phi(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts(it))); legend('$|\hat\phi|$');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor((RR),(ZZ),plt(si00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$r$'); ylabel('$z$'); title(sprintf('t=%.3d',Ts(it))); legend('$|S_i^{00}|$');
end
save_figure
end

%%
if 0
%% Show frame in kspace
tf = 20; [~,it] = min(abs(Ts-tf));
fig = figure; FIGNAME = ['krkz_frame',sprintf('_%.2d',JOBNUM)];
    subplot(221); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(PHI(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$|\hat\phi|$');
    subplot(222); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ni00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$|\hat n_i^{00}|$');
    subplot(223); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ne00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$|\hat n_e^{00}|$');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Si00(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_i^{00}$');
end
save_figure
end


if 0
%% Check time evolution of higher moments
p = 0; j = 0;
if strcmp(OUTPUTS.write_moments,'.true.') 

end
end
