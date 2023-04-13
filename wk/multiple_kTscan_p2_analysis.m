clrs_ = lines(10);
kN=2.22;

scantype = 'nuscan';
% scantype = 'kTscan';

switch scantype
    case 'kTscan'
%     resdir = 'paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/5x3x128x64x24_dp'; clr_ = clrs_(1,:); xname = '$\kappa_T (\kappa_N=2.22)$';
%     resdir = 'paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/7x4x128x64x24_dp'; clr_ = clrs_(2,:); xname = '$\kappa_T (\kappa_N=2.22)$';
    resdir = 'paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/9x5x128x64x24_dp'; clr_ = clrs_(3,:); xname = '$\kappa_T (\kappa_N=2.22)$';
    case 'nuscan'
%     resdir = 'paper_2_GYAC23/collision_study/nuDGGK_scan_kT_5.3/5x3x128x64x24_dp'; clr_ = clrs_(1,:); xname = '$\nu (\kappa_T=5.3,\kappa_N=2.22)$';
%     resdir = 'paper_2_GYAC23/collision_study/nuDGGK_scan_kT_5.3/9x5x128x64x24_dp'; clr_ = clrs_(3,:); xname = '$\nu (\kappa_T=5.3,\kappa_N=2.22)$';

%     resdir = 'paper_2_GYAC23/collision_study/nuSGGK_scan_kT_5.3/5x3x128x64x24_dp'; clr_ = clrs_(1,:); xname = '$\nu (\kappa_T=5.3,\kappa_N=2.22)$';
    resdir = 'paper_2_GYAC23/collision_study/nuSGGK_scan_kT_5.3/9x5x128x64x24_dp'; clr_ = clrs_(3,:); xname = '$\nu (\kappa_T=5.3,\kappa_N=2.22)$';
end

Njobs = 4;

x = 0*(1:Njobs);
Qx_avg  = 0*(1:Njobs);
Qx_std  = 0*(1:Njobs);
Chi_avg = 0*(1:Njobs);
Chi_std = 0*(1:Njobs);
figure
datadir = ['/misc/gyacomo23_outputs/',resdir,'/'];
for i = 1:Njobs+1
    J0 = i-1; J1 = i-1;

    Nseg = 5;

    data    = compile_results_low_mem(data,datadir,J0,J1);
    Trange  = data.Ts0D(end)*[0.3 1.0];
    %
    [~,it0] = min(abs(Trange(1)  -data.Ts0D)); 
    [~,it1] = min(abs(Trange(end)-data.Ts0D)); 
    %
    if 0
        Qx      = data.HFLUX_X(it0:it1);
        Qxa_    = 0*(1:Nseg);
        for n = 1:Nseg
           ntseg = floor((it1-it0)/n);
           for m = 1:n 
            Qxa_(n) = Qxa_(n) + mean(Qx((1:ntseg)+(m-1)*ntseg));
           end
           Qxa_(n) = Qxa_(n)/n;
        end
        Qx_avg(i) = mean(Qxa_);
        Qx_std(i) = std(Qxa_);
    else
        Qx_avg(i) = mean(data.HFLUX_X(it0:it1));
        Qx_std(i) =  std(data.HFLUX_X(it0:it1));
    end
    Chi_avg(i) = Qx_avg(i)./data.inputs.K_T/data.inputs.K_N;
    Chi_std(i) = Qx_std(i)./data.inputs.K_T/data.inputs.K_N;
    switch scantype
        case 'nuscan'
        x(i) = data.inputs.NU;
        case 'kTscan'
        x(i) = data.inputs.K_T;
    end
end

% plot;
errorbar(x,Chi_avg,Chi_std,'DisplayName',data.paramshort,'color',clr_); hold on;
switch  scantype
    case 'nuscan'
    Lin1999 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Lin_1999_fig2.txt');
    plot(Lin1999(:,1),Lin1999(:,2),'--ok','DisplayName','Lin1999');
    case 'kTscan'
    Dim2000 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Dimits_2000_fig3_full_no_GF.txt');
    plot(Dim2000(:,1),Dim2000(:,2),'ok','DisplayName','Dimits 2000');
end
% plot(ITG_threshold*[1 1],[0 20],'-.','DisplayName','$\kappa_T^{crit}$',...
%     'color',clr_);
ylabel('$\chi$');
xlabel(xname);
ylim([0,5]);
legend('show');
