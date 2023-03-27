
function [ fig, res ] = statistical_transport_averaging( data, options )
scale = data.scale;
Trange  = options.T;
[~,it0] = min(abs(Trange(1)-data.Ts0D)); 
[~,it1] = min(abs(Trange(end)-data.Ts0D)); 
gamma   = data.PGAMMA_RI(it0:it1)*scale;
Qx      = data.HFLUX_X(it0:it1)*scale;
dt_const = numel(unique(round(diff(data.Ts0D(it0:it1))*100)))==1;
% if ~dt_const
%     disp('DT not const on given interval');
% else
    
    Ntot = (it1-it0)+1;

    transp_seg_avg = 1:Ntot;
    transp_seg_std = 1:Ntot;

    for Np = 1:Ntot % Loop on the number of segments
        transp_seg_avg(Np) = mean(gamma(1:Np));
        transp_seg_std(Np) = std(gamma(1:Np));
    end

    time_seg = (data.Ts0D(it0:it1)); 
    
    fig = 0;
if options.NPLOTS > 0
    fig = figure;
    subplot(211)
    plot(time_seg,transp_seg_avg,'-'); hold on;
    xlabel('Averaging time'); ylabel('$\langle\Gamma_x\rangle_{\tau}$');
    legend(['$\Gamma_x^\infty=$',sprintf('%2.2e',transp_seg_avg(end))])
    title(sprintf('Transport averaging from t=%2.2f',data.Ts0D(it0)));
   
    for Np = 1:Ntot % Loop on the number of segments
        transp_seg_avg(Np) = mean(Qx(1:Np));
    end

    
    subplot(212)

    plot(time_seg,transp_seg_avg,'-'); hold on;
    xlabel('Averaging time'); ylabel('$\langle Q_x\rangle_{\tau}$');
    legend(['$Q_x^\infty=$',sprintf('%2.2e',transp_seg_avg(end))])
end   
res.time_seg = time_seg;
res.Qx_t     = transp_seg_avg;
res.Gx_avg = mean(gamma);
res.Gx_std = std (gamma);
disp(['G_x=',sprintf('%2.2e',res.Gx_avg),'+-',sprintf('%2.2e',res.Gx_std)]);
res.Qx_avg = mean(Qx);
res.Qx_std = std (Qx);
disp(['G_x=',sprintf('%2.2e',res.Qx_avg),'+-',sprintf('%2.2e',res.Qx_std)]);
end

