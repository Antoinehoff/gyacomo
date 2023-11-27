function [ fig, arrays, R, Z] = plot_metric( data, options )
fig = 0;
% names = {'Jacobian','gradxB','gradyB','gradzB','gradz_coeff',...
%          'gxx','gxy','gxz','gyy','gyz','gzz','hatB','hatR','hatZ'};
names = {'gxx','gxy','gxz','gyy','gyz','gzz',...
         'hatB', 'dBdx', 'dBdy', 'dBdz',...
         'Jacobian','hatR','hatZ','gradz_coeff'};
geo_arrays = zeros(data.grids.Nz,numel(names));

for i_ = 1:numel(names)
    namae = names{i_};
    geo_arrays(:,i_) = h5read(data.outfilenames{end},['/data/metric/',namae])';
end
gxx = geo_arrays(:,1);
gxy = geo_arrays(:,2);
gxz = geo_arrays(:,3);
gyy = geo_arrays(:,4);
gyz = geo_arrays(:,5);
% gzz = geo_arrays(:,6);
B   = geo_arrays(:,7);
dxB = geo_arrays(:,8);
dyB = geo_arrays(:,9);
dzB = geo_arrays(:,10);
Jac = geo_arrays(:,11);

try
    NPLOT = options.SHOW_FLUXSURF + options.SHOW_METRICS + options.SHOW_CURVOP;
catch
    NPLOT = 2;
    options.SHOW_FLUXSURF = 1;
    options.SHOW_METRICS  = 0;
    options.SHOW_CURVOP   = 1;
end
if NPLOT > 0
    fig = figure; 
    if options.SHOW_METRICS
    subplot(3,NPLOT,1*NPLOT)
        for i = 1:6
        plot(data.grids.z, geo_arrays(:,i),'DisplayName',names{i}); hold on;
        end
        xlim([min(data.grids.z),max(data.grids.z)]); legend('show'); title('GYACOMO geometry');

    subplot(3,NPLOT,2*NPLOT)
        for i = 7:10
        plot(data.grids.z, geo_arrays(:,i),'DisplayName',names{i}); hold on;
        end
        xlim([min(data.grids.z),max(data.grids.z)]); legend('show');

    subplot(3,NPLOT,3*NPLOT)
        for i = 11:14
        plot(data.grids.z, geo_arrays(:,i),'DisplayName',names{i}); hold on;
        end
        xlim([min(data.grids.z),max(data.grids.z)]); legend('show');
        legend('Interpreter','none')
    end
    if options.SHOW_FLUXSURF
        subplot(1,NPLOT,1)
        R = [geo_arrays(:,12);geo_arrays(1,12)];
        Z = [geo_arrays(:,13);geo_arrays(1,13)];
    plot(R./data.inputs.EPS,Z./data.inputs.EPS,'-k');
    xlabel('R [m]'); ylabel('Z [m]');
    axis tight
    axis equal
    end
    if options.SHOW_CURVOP
        G1 = gxx.*gyy - gxy.*gxy;
        G2 = gxx.*gyz - gxy.*gxz;
        G3 = gxy.*gyz - gyy.*gxz;
        Ckx0 = -(dyB./B + G2./G1.*dzB./B); 
        C0ky =  (dxB./B - G3./G1.*dzB./B); 
        % subplot(1,NPLOT,2);
        %     plot(data.grids.z,Ckx0,'DisplayName','$\mathcal C(k_x,k_y=0)/ik_x$'); hold on
        %     plot(data.grids.z,C0ky,'DisplayName','$\mathcal C(k_x=0,k_y)/ik_y$');
        %     xlabel('$z$'); legend("show");
        subplot(1,NPLOT,2);
            % z = [data.grids.z; data.grids.z(1)]+pi; 
            z = data.grids.z;
            % rho = 1+0.5*(Ckx0)./max(abs(Ckx0));
            rho = Ckx0;
            % rho = [rho; rho(1)];
            minr = min(rho); maxr = max(rho);
            plot(z,rho,'DisplayName','$\mathcal C(k_x,k_y=0)/ik_x$'); hold on
            % rho = 1+0.5*(C0ky)./max(abs(C0ky));
            rho = C0ky;
            % rho = [rho; rho(1)];
            minr = min(minr,min(rho)); maxr = max(maxr,max(rho));
            plot(z,rho,'DisplayName','$\mathcal C(k_x=0,k_y)/ik_y$');
            % rho = 1+0.5*(B.^(-1)./Jac)./max(abs(B.^(-1)./Jac));
            rho = B.^(-1)./Jac;
            % rho = [rho; rho(1)];
            minr = min(minr,min(rho)); maxr = max(maxr,max(rho));
            plot(z,rho,'DisplayName','$\hat B^{-1}/J_{xyz}$');
            % polarplot(z,zeros(size(rho)),'-k','DisplayName','$0$')
            legend("show",'Location','south','FontSize',14);
            % rlim(1.1*[minr maxr]);
    end
end
%outputs
arrays = squeeze(geo_arrays(:,:));
R = [geo_arrays(:,12);geo_arrays(1,12)];
Z = [geo_arrays(:,13);geo_arrays(1,13)];
try
    if ~options.SHOWFIG
        close
    end
catch
end
end

