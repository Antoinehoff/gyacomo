function plot_fa_gene(OPTIONS)
folder = OPTIONS.folder;
TIMES  = OPTIONS.times;
specie = OPTIONS.specie;
PLT_FCT= OPTIONS.PLT_FCT;

file = 'coord.dat.h5';
vp = h5read([folder,file],'/coord/vp'); nvp = numel(vp);
mu = h5read([folder,file],'/coord/mu'); nmu = numel(mu);
z  = h5read([folder,file],'/coord/z');
[XX,SS] = meshgrid(mu,vp);

switch OPTIONS.Z
    case 'avg'
        zcomp_name = ' z-avg';
        zcomp = @(x) squeeze(mean(x,1));
    otherwise
        zcomp_name = [' z=',sprintf('%2.2f',z(OPTIONS.Z))];
        zcomp = @(x) squeeze(x(OPTIONS.Z,:,:));
end

[~,iv0] = min(abs(vp));
[~,im0] = min(abs(mu));

file = 'vsp.dat.h5';
time  = h5read([folder,file],'/vsp/time');

fig = figure;

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
    switch specie
        case 'e'
        FFa    = squeeze(Gdata(1,:,:,2));
        FFa    = abs(FFa)./max(max(abs(FFa)));
    subplot(1,2,1)
        plot(vp,FFa(:,im0)); hold on;
    subplot(1,2,2)
        plot(mu,FFa(iv0,:)); hold on;
        case 'i'
        FFa   = squeeze(Gdata(1,:,:,1));    
        FFa    = abs(FFa)./max(max(abs(FFa)));
    end

    subplot(1,2,1)
        plot(vp,FFa(:,im0)); hold on;
        legend(specie)
        xlabel('$v_\parallel, (\mu=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 

    subplot(1,2,2)
        plot(mu,FFa(iv0,:)); hold on;
        legend(specie)
        xlabel('$\mu, (v_\parallel=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 
else
switch specie
    case 'e'
    name  = '$f_e(v_\parallel,\mu_p)$';
    FFa    = zcomp(squeeze(Gdata));
    FFa    = abs(FFa)./max(max(abs(FFa)));
    switch PLT_FCT
        case 'contour'
            contour(SS,XX,FFa);
        case 'contourf'
            pclr = contourf(SS,XX,FFa);
        case 'pcolor'
            pclr = pcolor(SS,XX,FFa); set(pclr, 'edgecolor','none'); shading interp
    end
    case 'i'
    name  = '$f_i(v_\parallel,\mu_p)$';
    FFa    = zcomp(squeeze(Gdata));
    FFa    = abs(FFa)./max(max(abs(FFa)));
    switch PLT_FCT
        case 'contour'
            contour(SS,XX,FFa,(nvp+nmu)/2);
        case 'contourf'
            pclr = contourf(SS,XX,FFa,(nvp+nmu)/2);
        case 'pcolor'
            pclr = pcolor(SS,XX,FFa); set(pclr, 'edgecolor','none'); shading interp
    end
end
if numel(TIMES) == 1
    title(['Gene ',name,zcomp_name,', $t = ',sprintf('%2.2f',time(it)),'$']);
else
    title(['Gene ',name,zcomp_name,', averaged $t\in$[',sprintf('%2.2f',G_t(1)),',',sprintf('%2.2f',G_t(end)),']']);
end

clear FFi FFe tmp Gdata
end