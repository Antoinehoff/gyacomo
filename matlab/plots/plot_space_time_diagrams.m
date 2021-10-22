[~,itstart] = min(abs(Ts3D-tstart));
[~,itend]   = min(abs(Ts3D-tend));
trange = itstart:itend;
[TY,TX] = meshgrid(x,Ts3D(trange));
fig = figure; FIGNAME = ['space_time','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(211)
        pclr = pcolor(TX,TY,squeeze(mean(dens_i(:,:,1,trange).*dyphi(:,:,1,trange),2))'); set(pclr, 'edgecolor','none'); colorbar;
%         pclr = pcolor(TX,TY,squeeze(mean(ni00(:,:,1,trange).*dyphi(:,:,1,trange),2))'); set(pclr, 'edgecolor','none',...
%             'DisplayName','$\langle n_i\partial_z\phi\rangle_z$'); colorbar;
        shading interp
        colormap hot;
%         caxis([0.0,0.05*max(max(mean(ni00(:,:,its2D:ite2D).*dyphi(:,:,1,its2D:ite2D),2)))]);
        caxis([0.0,cmax]); c = colorbar; c.Label.String ='\langle\Gamma_{x}\rangle_{z}';
         xticks([]); ylabel('$x/\rho_s$')
%         legend('Radial part. transport $\langle n_i\partial_z\phi\rangle_z$')
        title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\kappa_N=$',num2str(K_N),...
        ', $L=',num2str(L),'$, $N=',num2str(Nx),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
    subplot(212)
        pclr = pcolor(TX,TY,squeeze(mean(dxphi(:,:,1,trange),2))'); set(pclr, 'edgecolor','none',...
            'DisplayName','$\langle \partial_r\phi\rangle_z$'); colorbar;
        fieldmax = max(max(mean(abs(dxphi(:,:,1,its2D:ite2D)),2)));
        caxis([-fieldmax,fieldmax]); c = colorbar; c.Label.String ='\langle v_{E\times B,z}\rangle_z';
        colormap(bluewhitered)
        xlabel('$t c_s/R$'), ylabel('$x/\rho_s$')
%         legend('Zonal flow $\langle \partial_r\phi\rangle_z$')
save_figure
