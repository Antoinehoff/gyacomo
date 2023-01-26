function [FIGURE] =  plot_fa_gene(OPTIONS)
folder = OPTIONS.folder;
TIMES  = OPTIONS.T;
SPECIES = OPTIONS.SPECIES;
PLT_FCT= OPTIONS.PLT_FCT;

file = 'coord.dat.h5';
vp = h5read([folder,file],'/coord/vp'); nvp = numel(vp);
mu = h5read([folder,file],'/coord/mu'); nmu = numel(mu);
z  = h5read([folder,file],'/coord/z');
[XX,SS] = meshgrid(mu,vp);

switch OPTIONS.iz
    case 'avg'
        zcomp_name = ' z-avg';
        zcomp = @(x) squeeze(mean(x,1));
    otherwise
        zcomp_name = [' z=',sprintf('%2.2f',z(OPTIONS.iz))];
        zcomp = @(x) squeeze(x(OPTIONS.iz,:,:));
end

[~,iv0] = min(abs(vp));
[~,im0] = min(abs(mu));

file = 'vsp.dat.h5';
time  = h5read([folder,file],'/vsp/time');

FIGURE.fig = figure;

G_t = [];
Gdata = 0;
for T = TIMES
[~, it] = min(abs(time-T));
tmp   = h5read([folder,file],['/vsp/',OPTIONS.FIELD,'/',sprintf('%10.10d',it-1)]);
Gdata = Gdata + tmp;
G_t = [G_t time(it)];
end
Gdata = Gdata ./ numel(TIMES);

if OPTIONS.ONED
    switch SPECIES
        case 'e'
        FFa    = zcomp(Gdata(:,:,:,2));
        FFa    = abs(FFa)./max(max(abs(FFa)));
    FIGURE.ax1 = subplot(1,2,1,'parent',FIGURE.fig);
        plot(vp,FFa(:,im0)); hold on;
    FIGURE.ax2 = subplot(1,2,2,'parent',FIGURE.fig);
        plot(mu,FFa(iv0,:)); hold on;
        case 'i'
        FFa   = zcomp(Gdata(:,:,:,1));    
        FFa    = abs(FFa)./max(max(abs(FFa)));
    end

    FIGURE.ax1 = subplot(1,2,1,'parent',FIGURE.fig);
        plot(vp,FFa(:,im0)); hold on;
        legend(SPECIES)
        xlabel('$v_\parallel, (\mu=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 

    FIGURE.ax2 = subplot(1,2,2,'parent',FIGURE.fig);
        plot(mu,FFa(iv0,:)); hold on;
        legend(SPECIES)
        xlabel('$\mu, (v_\parallel=0)$'); %ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 
        xticks(FIGURE.ax2,[]);
else
    FIGURE.ax1 = subplot(1,1,1,'parent',FIGURE.fig);
switch SPECIES
    case 'e'
    name  = '$f_e(v_\parallel,\mu_p)$';
    FFa    = zcomp(squeeze(Gdata(:,:,:,2)));
    FFa    = abs(FFa)./max(max(abs(FFa)));
    switch PLT_FCT
        case 'contour'
            contour(SS,XX,FFa);
        case 'contourf'
            pclr = contourf(SS,XX,FFa);
        case 'pcolor'
            pclr = pcolor(SS,XX,FFa); set(pclr, 'edgecolor','none'); shading interp
        case 'surf'
            surf(SS,XX,FFa);
        case 'surfvv'
            surf([SS; SS(:,end-1:-1:1)],sqrt([XX; XX(:,end-1:-1:1)]),[FFa; FFa(:,end-1:-1:1)]);
    end
    case 'i'
    name  = '$f_i(v_\parallel,\mu_p)$';
    FFa    = zcomp(squeeze(Gdata(:,:,:,1)));
    FFa    = abs(FFa)./max(max(abs(FFa)));
    switch PLT_FCT
        case 'contour'
            contour(SS,XX,FFa,(nvp+nmu)/2);
        case 'contourf'
            pclr = contourf(SS,XX,FFa,(nvp+nmu)/2);
        case 'pcolor'
            pclr = pcolor(SS,XX,FFa); set(pclr, 'edgecolor','none'); shading interp
        case 'surf'
            surf(SS,XX,FFa);
        case 'surfvv'
            surf([SS(:,end:-1:1) SS ],[-sqrt(XX(:,end:-1:1)) sqrt(XX)],[FFa(:,end:-1:1) FFa]);
            xlabel('$v_\parallel$'); ylabel('$v_\perp$');
            xlim([min(vp) max(vp)]);
            ylim(sqrt(max(mu))*[-1 1]);
    end
end
if numel(TIMES) == 1
    title(['Gene ',name,zcomp_name,', $t = ',sprintf('%2.2f',time(it)),'$']);
else
    title(['Gene ',name,zcomp_name,', averaged $t\in$[',sprintf('%2.2f',G_t(1)),',',sprintf('%2.2f',G_t(end)),']']);
end

clear FFi FFe tmp Gdata
end