% Get the current directory (wk)
curdir  = pwd;
NCONTOUR = 0;
% give ref value and param names
REFVAL= 0;
% normalize to the max all data
NORM_ALL= 0;
% normalize to the max each line
NORM_LIN= 0;
% normalize to the max each column
NORM_COL= 0;
% Get and plot the fluxsurface
GETFLUXSURFACE = 0;

% partition= '../results/paper_3/';
% Get the scan directory
switch 6
    case 1 % delta kappa scan
        casename = 'DTT rho85';
        partition= '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/';
        scandir = 'P2_J1_delta_kappa_scan'; scanname= '(2,1)';
        % scandir = 'P4_J2_delta_kappa_scan'; scanname= '(4,2)';
        nml1 = 'GEOMETRY'; pnam1 = '$\delta$'; attr1 = 'delta'; pref1 = 0.23; scale1 =1.0;
        nml2 = 'GEOMETRY'; pnam2 = '$\kappa$'; attr2 = 'kappa'; pref2 = 1.53; scale2 =1.0;
        t1 = 50; t2 = 150;
    case 2 % shear safety factor scan
        casename = 'DTT rho85';
        partition= '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/';
        scandir = 'P2_J1_PT_sfact_shear_scan'; scanname= '(2,1)';
        % scandir = 'P2_J1_NT_sfact_shear_scan'; scanname= '(2,1)';
        nml1 = 'GEOMETRY'; pnam1 = '$\hat s$'; attr1 = 'shear'; pref1 = 3.63; scale1 =1.0;
        nml2 = 'GEOMETRY'; pnam2 = '$q_0$';    attr2 = 'q0';    pref2 =-2.15; scale2 =1.0;
        t1 = 50; t2 = 150;
    case 3
        casename = 'DTT rho85';
        partition= '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/';
        % scandir = 'P2_J1_delta_nuDGGK_scan';       scanname= 'DG (2,1)';
        % scandir = 'P4_J2_delta_nuDGGK_scan';      scanname= 'DG (4,2)';
        % scandir = 'P8_J4_delta_nuDGGK_conv_test'; scanname= 'DG (8,4)';
        % scandir = 'P2_J1_delta_nuSGGK_scan';      scanname= 'SG (2,1)';
        % scandir = 'P4_J2_delta_nuSGGK_scan';      scanname= 'SG (4,2)';
        % scandir = 'P8_J4_delta_nuSGGK_conv_test'; scanname= 'SG (8,4)';
        % scandir = 'P4_J2_delta_nuSGGKii_scan';    scanname= 'SGii (4,2)';
        scandir = 'P2_J1_delta_nuLDGK_scan';      scanname= 'LD (2,1)';
        % scandir = 'P4_J2_delta_nuLDGK_scan';      scanname= 'LD (4,2)';
        nml1 = 'GEOMETRY'; pnam1 = '$\delta$'; attr1 = 'delta'; pref1 = 0.23; scale1 =1.0;
        nml2 = 'MODEL';    pnam2 = '$\nu$';    attr2 = 'nu';    pref2 = 0.5;  scale2 =1.0;
        t1 = 50; t2 = 150;
    case 4 % shear delta scan
        casename = 'DIIID rho95';
        partition= '/misc/gyacomo23_outputs/paper_3/DIIID_cold_ions_rho95_geom_scan/';
        scandir = '3x2x192x48x32_RT_2500_delta_shear_scan'; scanname= 'CI DG RT 2500';
        % scandir = '3x2x256x64x48_delta_shear_scan';
        nml1 = 'GEOMETRY'; pnam1 = '$\delta$'; attr1 = 'delta'; pref1 = 0;   scale1 =1.0;
        nml2 = 'GEOMETRY'; pnam2 = '$\hat s$'; attr2 = 'shear'; pref2 = 0.8; scale2 =1.0;
        t1 = 50; t2 = 150;
    case 5 % delta K_T tau=1
        casename = 'DIIID rho95 $\tau=1$';
        partition= '../results/paper_3/DIIID_tau_1_rho95_geom_scan/';  
        % scandir = '3x2x192x48x32_delta_RT_scan'; scanname= '(2,1)';
        scandir = '2_1_delta_RT_scan'; scanname= '(2,1)';
        % scandir = '5x3x192x48x32_delta_RT_scan'; scanname= '(4,2)';
        % scandir = 'delta_RT_scan_PJ_21'; scanname= '(2,1)';
        nml1 = 'GEOMETRY'; pnam1 = '$\delta$'; attr1 = 'delta'; pref1 = 0; scale1 =1.0;
        nml2 = 'SPECIES'; pnam2 = '$R_0/L_T\times T_i/T_e$'; attr2 = 'K_T_'; pref2 = 0.8; scale2 =1.0;
        t1 = 200; t2 = 500;
    case 6 % delta K_T cold ions
        casename = 'DIIID rho95 $\tau=10^{-3}$';
        partition= '../results/paper_3/DIIID_cold_ions_rho95_geom_scan/'; 
        scandir = '3x2x192x48x32_delta_RT_scan'; scanname= '(2,1)';
        nml1 = 'GEOMETRY'; pnam1 = '$\delta$'; attr1 = 'delta'; pref1 = 0; scale1 =1.0;
        nml2 = 'SPECIES'; pnam2 = '$R_0/L_T\times T_i/T_e$'; attr2 = 'K_T_'; pref2 = 0.8; scale2 =1e3/2;
        t1 = 150; t2 = 280;
   case 7 % delta s_delta
        casename = 'DIIID rho95 $\tau=10^{-3}$';
        partition= '/misc/gyacomo23_outputs/paper_3/DIIID_cold_ions_rho95_geom_scan/';
        scandir  = '3x2x192x48x32_RT_1000_delta_sdelta_scan'; scanname= 'RT=1000 (2,1)';
        % scandir = '';
        nml1 = 'GEOMETRY'; pnam1 = '$\delta$'; attr1 = 'delta'; pref1 = 0; scale1 =1.0;
        nml2 = 'GEOMETRY'; pnam2 = '$s_\delta$'; attr2 = 's_delta';  pref2 = 0.8; scale2 =1.0;
        t1 = 200; t2 = 295;
    case 8 % eps q0
        casename = 'DIIID rho95 $\tau=10^{-3}$';
        partition= '/misc/gyacomo23_outputs/paper_3/DIIID_cold_ions_rho95_geom_scan/';
        scandir  = '3x2x192x48x32_RT_1000_eps_q0_scan/PT'; scanname= 'PT, RT=1000 (2,1)';
        % scandir  = '3x2x192x48x32_RT_1000_eps_q0_scan/NT'; scanname= 'NT, RT=1000 (2,1)';
        % scandir = '';
        nml1 = 'GEOMETRY'; pnam1 = '$\epsilon$'; attr1 = 'eps'; pref1 = 0; scale1 =1.0;
        nml2 = 'GEOMETRY'; pnam2 = '$q_0$';      attr2 = 'q0';      pref2 = 0; scale2 =1.0;
        t1 = 200; t2 = 400;
    case 9 % CBC Dimits shift
        casename = 'HEL CBC';
        partition= '/misc/gyacomo23_outputs/paper_3/HEL_CBC/';
        scandir  = '128x32x24'; scanname= 'CBC HEL';
        nml1 = 'SPECIES'; pnam1 = '$R_T$'; attr1 = 'k_T_'; pref1 = 0; scale1 =1.0;
        nml2 = 'SPECIES'; pnam2 = '$R_N$'; attr2 = 'k_N_'; pref2 = 0; scale2 =1.0;
        t1 = 1000; t2 = 2000;
end 
scanname= [casename scanname];
scandir = [partition,scandir,'/']; 
% Get a list of all items in the current directory
contents = dir(scandir);

% Iterate through the contents
Qxavg = []; Qxerr = []; para1 = []; para2 = []; R = []; Z = [];
Qxt = struct();
for i = 1:length(contents)
    % Check if the item is a directory and not '.' or '..'
    if contents(i).isdir && ~strcmp(contents(i).name, '.') ...
            && ~strcmp(contents(i).name, '..')
        % Get and display the name of the subdirectory
        subdir = [scandir,contents(i).name];
        disp(['Subdirectory: ' contents(i).name]);
        % Get parameters
        param = read_namelist([subdir,'/fort_00.90']);
        para1 = [para1 param.(nml1).(attr1)];
        para2 = [para2 param.(nml2).(attr2)];        
        % Now you are in the subdirectory. You can perform operations here.
        [t_all, Pxi_all, Qxi_all, Pxe_all, Qxe_all] = read_flux_out_XX(subdir);
        if(numel(Qxe_all) > 1)
            Qxtot = Qxi_all+Qxe_all;
        else
            Qxtot = Qxi_all;
        end
        Qxt.(['dat_',num2str(i)])      = struct();
        Qxt.(['dat_',num2str(i)]).Qx   = Qxtot;
        Qxt.(['dat_',num2str(i)]).t    = t_all;
        Qxt.(['dat_',num2str(i)]).name = contents(i).name;
        if(numel(t_all) > 1)
          disp(num2str(t_all(end)))
            [~,it1]  = min(abs(t_all-t1));
            [~,it2]  = min(abs(t_all-t2));
            steady_slice = it1:it2;
            if(t_all(end) >= t2)
                [fullAvg,sliceAvg,sliceErr] = sliceAverage(Qxtot(steady_slice),3);
                Qxavg = [Qxavg fullAvg];
                Qxerr = [Qxerr mean(sliceErr)];
            else
                Qxavg = [Qxavg nan];
                Qxerr = [Qxerr nan];
            end
        else
                            Qxavg = [Qxavg nan];
                Qxerr = [Qxerr nan];
        end
    end
    if GETFLUXSURFACE
        data = load([subdir,'/RZ.txt']);
        R_ = data(:, 1);
        Z_ = data(:, 2);
        R_ = [R_;R_(1)]'; Z_  = [Z_;Z_(1)]';
        R  = [R ; R_];    Z   = [Z ; Z_];
    end

end
if 0
%% plot time traces
attr = fieldnames(Qxt);
Nsim = numel(attr);
figure
% compute growth at the begining
tw = [10 40];
gr = 1:Nsim; err = 1:Nsim;
for i = 1:1:Nsim
    tmp_ = Qxt.(attr{i});
    t = tmp_.t;
    y = tmp_.Qx;  
    plot(t,y,'DisplayName',tmp_.name); hold on;
    [~,it1] = min(abs(t-tw(1)));
    [~,it2] = min(abs(t-tw(2)));
    [gr_, err_] = compute_growth(t(it1:it2),y(it1:it2));
    gr(i) = gr_; err(i) = err_;
end
%%
toplot = real(reshape(gr,sz))';
toplot = toplot(idx1,idx2);

figure
imagesc_custom(xx_,yy_,toplot); hold on
end
%% reshaping, sorting and plotting
p1 = unique(para1)/scale1;
p2 = unique(para2)/scale2;
N1 = numel(p1);
N2 = numel(p2);

if para1(1) == para1(2)
    sz = [N2 N1];
    TRANSPOSE = 1;
else
    sz = [N1 N2];
    TRANSPOSE = 0;
end

Zavg = reshape(Qxavg,sz);
Zerr = reshape(Qxerr,sz);
XX   = reshape(para1/scale1,sz);
YY   = reshape(para2/scale2,sz);

if TRANSPOSE
    Zavg = Zavg';
    Zerr = Zerr';
    XX = XX';
    YY = YY';
end

[~,idx1] = sort(XX(:,1));
[~,idx2] = sort(YY(1,:));
Zavg = Zavg(idx1,idx2);
Zerr = Zerr(idx1,idx2);
XX   = XX(idx1,idx2);
YY   = YY(idx1,idx2);

% compute the 
if REFVAL
    if NORM_ALL
    Qxname = '$\bar Q_{tot}/\bar Q_{max}[\%]$';
        [tmp,iref1] = max(Zavg);
        [~,  iref2] = max(tmp);
        iref1 = iref1(iref2);
    else
    Qxname = '$\langle (Q_{tot}-Q_{ref})/Q_{ref} \rangle_t[\%]$';
        if pref1 ~= 999
            [~,iref1] = min(abs(XX(:,1)-pref1));
        else
            iref1     = 1:N1;
        end
        if pref2 ~= 999
            [~,iref2] = min(abs(YY(1,:)-pref2));
        else
            iref2     = 1:N2;
        end
    end
    iref1     = ones(N1,1).*iref1;
    iref2     = ones(N2,1).*iref2;
    xref  = XX(iref1,iref2);
    yref  = YY(iref1,iref2);
    Qxref = Zavg(iref1,iref2);
    Qrefname = ['$Q_{ref}=$',num2str(Qxref(1,1))];
else
    Qref = 1;
    if NORM_LIN
        Qxname = '$\bar Q_{tot}/\bar Q_{max}[\%]$, per line';
        for il = 1:sz(1)
            maxline = max(Zavg(:,il));
            Zavg(:,il) = Zavg(:,il)./maxline;
            Zerr(:,il) = Zerr(:,il)./maxline;
        end
    elseif NORM_COL
        Qxname = '$\bar Q_{tot}/\bar Q_{max}[\%]$, per column';
        for ic = 1:sz(2)
            maxcol = max(Zavg(ic,:));
            Zavg(ic,:) = Zavg(ic,:)./maxcol;
            Zerr(ic,:) = Zerr(ic,:)./maxcol;
        end
    else
        Qxname = '$\langle Q_{tot} \rangle_t$';
    end
end

% Figure
figure
subplot(1,2,1)
[xx_,yy_] = meshgrid(XX(:,1),YY(1,:));
if REFVAL
    if NORM_ALL || NORM_COL || NORM_LIN
        toplot = (Zavg./Qxref)'
        CLIM = [0 1];
    else
        toplot = ((Zavg-Qxref)./Qxref * 100)';
        CLIM = 'auto';
    end
else
    toplot = Zavg';
    CLIM = 'auto';
end
if NCONTOUR <= 0
    imagesc_custom(xx_,yy_,toplot); hold on
else
    contourf(XX(:,1),YY(1,:),Zavg',NCONTOUR); hold on
end
if REFVAL && ~((pref1==999) || (pref2==999))
    plot(xref(1,1),yref(1,1),'xk','MarkerSize',14,'DisplayName',Qrefname)
    legend('show')
end
xlabel(pnam1); ylabel(pnam2);
title(scanname)
colormap(bluewhitered); colorbar; clim(CLIM);
if ~REFVAL
    colormap(jet); 
end
subplot(1,2,2)
clrs = jet(N2);
for i = 1:N2
    errorbar(XX(:,i),Zavg(:,i),Zerr(:,i),...
        'DisplayName',[pnam2,'=',num2str(p2(i))],...
        'Color',clrs(i,:));
    hold on;
end
if REFVAL && ~((pref1==999) || (pref2==999))
    plot(xref(1,1),0,'xk','MarkerSize',14,'DisplayName',Qrefname)
end
grid on
xlabel(pnam1); ylabel('$\langle Q_{tot} \rangle_t$');
legend('show','Location','northwest');
title([param.COLLISION.collision_model{1}, ...
    ', $(P,J)=(',num2str(param.GRID.pmax),',',num2str(param.GRID.jmax),')$'])
