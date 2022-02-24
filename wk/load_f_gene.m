% folder = '/misc/gene_results/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-2_mu_1e-2_SGDK_36x20/';
folder = '/misc/gene_results/HP_fig_2a_mu_1e-2/';
% folder = '/misc/gene_results/HP_fig_2b_mu_5e-2/';
% folder = '/misc/gene_results/HP_fig_2b_mu_1e-1/';
% folder = '/misc/gene_results/HP_fig_2c_mu_5e-2/';
% folder = '/misc/gene_results/HP_fig_2c_mu_1e-2_muv_1e-1/';
% folder = '/misc/gene_results/HP_fig_2c_gyroLES/';
% folder = '/misc/gene_results/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-2_mu_1e-2_SGDK_128x32x48x32/';
% coordData=loadCoord(folder,-1,1,0,'ions');
% info=dat_file_parser(folder,-1);
% 
% chkpt_file=[folder,'vsp.dat'];
% fid=fopen(chkpt_file);
% 
% [~,~,gene_data]=GetFortranStep(fid,'DOUBLE',8,'LITTLE',1);
% fclose(fid);
% 
% nx = info.box.nx0;
% ny = info.box.nky0;
% nz = info.box.nz0;
% nv = info.box.nv0;
% nw = info.box.nw0;
% ns = numel(info.species);
% Lw = info.box.lw;
% Lv = info.box.lv;
% s    = linspace(-Lv,Lv,nv);
% x    = linspace(0,Lw,nw);
% [XX,SS] = meshgrid(x,s);
% 
% gridsize=(nv*nw*nz);
% Gdata =reshape(gene_data((4*gridsize+1):(5*gridsize)),[nz nv nw]);

file = 'coord.dat.h5';
vp = h5read([folder,file],'/coord/vp');
mu = h5read([folder,file],'/coord/mu');
[XX,SS] = meshgrid(mu,vp);

[~,iv0] = min(abs(vp));
[~,im0] = min(abs(mu));

file = 'vsp.dat.h5';
time  = h5read([folder,file],'/vsp/time');

TIMES = 500:1000;

fig = figure;

G_t = [];
Gdata = 0;
for T = TIMES
[~, it] = min(abs(time-T));
tmp   = h5read([folder,file],['/vsp/<f_>/',sprintf('%10.10d',it-1)]);
Gdata = Gdata + tmp;
G_t = [G_t time(it)];
end
Gdata = Gdata ./ numel(TIMES);

if 0
    FFe    = squeeze(Gdata(1,:,:,2));
    FFe    = abs(FFe)./max(max(abs(FFe)));
subplot(1,2,1)
    plot(vp,FFe(:,im0)); hold on;
subplot(1,2,2)
    plot(mu,FFe(iv0,:)); hold on;

    FFi   = squeeze(Gdata(1,:,:,1));    
    FFi    = abs(FFi)./max(max(abs(FFi)));
subplot(1,2,1)
    plot(vp,FFi(:,im0)); hold on;
    legend('e','i')
    xlabel('$v_\parallel, (\mu=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 
    
subplot(1,2,2)
    plot(mu,FFi(iv0,:)); hold on;
    legend('e','i')
    xlabel('$\mu, (v_\parallel=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 

else
subplot(1,2,1)
    FFe    = squeeze(Gdata(1,:,:,1));
    FFe    = abs(FFe)./max(max(abs(FFe)));
    contour(SS,XX,FFe,128);
%     pclr = pcolor(SS,XX,FFe); set(pclr, 'edgecolor','none'); shading interp
subplot(1,2,2)
    FFi   = squeeze(Gdata(1,:,:,2));    
    FFi    = abs(FFi)./max(max(abs(FFi)));
    contour(SS,XX,FFi,128);
%     pclr = pcolor(SS,XX,FFi); set(pclr, 'edgecolor','none'); shading interp
end

subplot(1,2,1)
if numel(TIMES) == 1
    title(['Gene, $t = ',num2str(time(it)),'$']);
else
    title([' Gene, average $t\in$[',num2str(G_t(1)),',',num2str(G_t(end)),']']);
end

clear FFi FFe tmp Gdata
    