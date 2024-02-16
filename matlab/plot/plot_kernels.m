if 0
%% Kernels
kmax=5;
nmax=6;
kp = linspace(0,kmax,100);


figure
for n_ = 0:nmax
    plot(kp,kernel(n_,kp),'DisplayName',['$\mathcal{K}_{',num2str(n_),'}$']);hold on;
end
ylim_ = ylim;
plot(kp(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
% J0 = besselj(0,kx);
% plot(kx,J0,'-r','DisplayName','$J_0$');
legend('show'); xlabel('k'); title('Kernel functions $\mathcal{K}_n$');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
%% Bessels and approx
wperp = linspace(0,1.5,4);
nmax1=5;
nmax2=10;
kmax=7;
figure
for id2 = 1:4
subplot(2,2,id2)
    v_ = wperp(id2);
    kp = linspace(0,kmax,100);

    J0 = besselj(0,kp*v_);
    A1 = 1 - kp.^2*v_^2/4;
    K1 = zeros(size(kp));
    K2 = zeros(size(kp));
    for n_ = 0:nmax1
        K1 = K1 + kernel(n_,kp).*polyval(LaguerrePoly(n_),v_^2);
    end
    for n_ = 0:nmax2
        K2 = K2 + kernel(n_,kp).*polyval(LaguerrePoly(n_),v_^2);
    end
    plot(kp,J0,'-k','DisplayName','$J_0(kv)$'); hold on;
    plot(kp,A1,'-r','DisplayName','$1 - k^2 v^2/4$');
    plot(kp,K1,'--b','DisplayName',['$\sum_{n=0}^{',num2str(nmax1),'}\mathcal K_n(k) L^n(v)$']);
    plot(kp,K2,'-b','DisplayName',['$\sum_{n=0}^{',num2str(nmax2),'}\mathcal K_n(k) L^n(v)$']);
    ylim_ = [-0.5, 1.0];
    plot(kp(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
    ylim(ylim_); xlabel('$k$')
    legend('show'); grid on; title(['$v = ',num2str(v_),'$'])
end
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%% from Sum_n Kn x Ln x Lj to Sum_n Kn Sum_s dnjs x Ls
wperp = [0 0.2 0.4];
kp    = linspace(0,3,50);
Jmax  = 8;
js     = 4;
fig = figure; set(gcf, 'Position',  [100, 100, numel(wperp)*300, 250])
%     suptitle(['$J_{ma}=',num2str(Jmax),', j=',num2str(j),'$'])
Kn = @(x__,y__) kernel(x__,y__);
% dnjs_ = @(n,j,s) dnjs(n,j,s);
dnjs_ = @(n,j,s) dnjs_GYAC(n+1,j+1,s+1);
for id2 = 1:numel(wperp)
subplot(1,numel(wperp),id2)
    v_ = wperp(id2);
    % J0 x Lj
    J0Lj_ = besselj(0,kp*v_)*polyval(LaguerrePoly(js),v_^2);
    % Taylor exp of J0
    A1 = (1 - kp.^2*v_^2/4)*polyval(LaguerrePoly(js),v_^2);
    % Sum_n^Nmax Kn x Ln x Lj
    KLL = zeros(size(kp));
    for n_ = 0:Jmax
        KLL = KLL + Kn(n_,kp).*polyval(LaguerrePoly(n_),v_^2);
    end
    KLL = KLL.*polyval(LaguerrePoly(js),v_^2);
    % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
    KdL1_ = zeros(size(kp));
    for n_ = 0:Jmax-js
        tmp_ = 0;
        for s_ = 0:n_+js
            tmp_ = tmp_ + dnjs_(n_,js,s_)*polyval(LaguerrePoly(s_),v_^2);
        end
        KdL1_ = KdL1_ + Kn(n_,kp).*tmp_;
    end
    % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
    KdL2_ = zeros(size(kp));
    for n_ = 0:Jmax
        tmp_ = 0;
        for s_ = 0:min(Jmax,n_+js)
            tmp_ = tmp_ + dnjs_(n_,js,s_)*polyval(LaguerrePoly(s_),v_^2);
        end
        KdL2_ = KdL2_ + Kn(n_,kp).*tmp_;
    end
    
    % plots
    plot(kp,J0Lj_,'-k','DisplayName','$J_0 L_j$');
    ylim_ = ylim;
    plot(kp,KLL,'-','DisplayName',['$\sum_{n=0}^{J}\mathcal K_n L_n L_j$']); hold on;
    plot(kp,KdL2_,':o','DisplayName',['$\sum_{n=0}^{J}\mathcal K_n \sum_{s=0}^{\min(J,n+j)} d_{njs} L_s$']);
    plot(kp,J0Lj_,'-k','DisplayName','$J_0 L_j$'); hold on;
    ylim_ = ylim;
    plot(kp,A1,':s','DisplayName','$(1 - l_\perp) L_j$'); hold on;
    plot(kp,KdL1_,':x','MarkerSize',10,...
        'DisplayName',['$\sum_{n=0}^{J-j}\mathcal K_n \sum_{s=0}^{n+j} d_{njs} L_s$']);
    plot(kp,J0Lj_,'-k','DisplayName','$J_0 L_j$'); hold on;
%     plot(kx(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
    ylim(ylim_); 
    xlabel('$k_{\perp}\rho_s$')
    % legend('show'); 
    grid on; title(['$w_\perp = ',num2str(v_),'$',' $J_{max}=',num2str(Jmax),', j=',num2str(js),'$'])
end
FIGNAME = ['/home/ahoffman/Dropbox/Applications/Overleaf/Hermite-Laguerre Z-pinch (HeLaZ Doc.)/figures/J0Lj_approx_J','_',num2str(Jmax),'_j_',num2str(js),'.eps'];
% saveas(fig,FIGNAME,'epsc');
end
if 0
%% Test Lj Ln = sum_s dnjs Ls
js = 4;
n = 4;
GYAC = 1;
smax = n+js;
wperp = linspace(0,5,50);
LjLn        = zeros(size(wperp));
dnjsLs      = zeros(size(wperp));
dnjsLs_bin  = zeros(size(wperp));
dnjsLs_stir = zeros(size(wperp));
if GYAC
    dnjsLs_GYAC = zeros(size(wperp));
    % Specify the filename
    filename = '/home/ahoffman/gyacomo/results/dbg/output.txt';
    % Use the load function to read data from the text file
    data = load(filename);
    % Extract columns
    indices = data(:, 1:3);
    values = data(:, 4);
    is_ = unique(indices(:,1));
    N_  = numel(is_);
    dnjs_GYAC = zeros(N_,N_,N_);
    for id2 = 1:numel(values)
        in_ = indices(id2,1);
        ij_ = indices(id2,2);
        is_ = indices(id2,3);
        dnjs_GYAC(in_,ij_,is_) = values(id2);
    end
end
  
id2=1;
for x_ = wperp
    LjLn(id2) = polyval(LaguerrePoly(js),x_)*polyval(LaguerrePoly(n),x_);
    dnjsLs(id2)      = 0;
    dnjsLs_bin(id2)  = 0;
    dnjsLs_stir(id2) = 0;
    for s_ = 0:smax
        dnjsLs(id2)      = dnjsLs(id2)      + dnjs(n,js,s_)*polyval(LaguerrePoly(s_),x_);
        dnjsLs_bin(id2)  = dnjsLs_bin(id2)  + dnjs_fact(n,js,s_)*polyval(LaguerrePoly(s_),x_);
        dnjsLs_stir(id2) = dnjsLs_stir(id2) + dnjs_stir(n,js,s_)*polyval(LaguerrePoly(s_),x_);
        if GYAC
        dnjsLs_GYAC(id2) = dnjsLs_GYAC(id2) + dnjs_GYAC(n+1,js+1,s_+1)*polyval(LaguerrePoly(s_),x_);
        end
    end
    id2 = id2+1;
end

figure
plot(wperp,LjLn,'DisplayName',['$L_',num2str(js),' L_',num2str(n),'$']); hold on;
plot(wperp,dnjsLs,'d','DisplayName',...
    ['$\sum_{s=0}^{',num2str(smax),'} d^{coeff}_{',num2str(n),',',num2str(js),',s}L_s$'])
plot(wperp,dnjsLs_bin,'o','DisplayName','$\sum_s d^{fact}_{njs}L_s$')
plot(wperp,dnjsLs_stir/dnjsLs_stir(1),'x','DisplayName','$\sum_s d^{stir}_{njs}L_s$')
plot(wperp,dnjsLs_GYAC,'x','DisplayName','$GYAC$')
xlabel('$w_\perp$');
% We see that starting arround j = 18, n = 0, the stirling formula is the only
% approximation that does not diverge from the target function. However it
% performs badly for non 0 n...
end

if 0
%%
    wperp = linspace(0,3,128);
    kp    = linspace(0,3,128);
    Jmax  = 8;
    js     = Jmax/2;[0 2:Jmax];
    err_KLL  = 1:numel(js);
    err_KdL1 = 1:numel(js);
    err_KdL2 = 1:numel(js);
    err_T1   = 1:numel(js);
    J0Lj     = zeros(numel(kp),numel(wperp));
    KdL1     = zeros(numel(kp),numel(wperp));
    KdL2     = zeros(numel(kp),numel(wperp));
    T1       = zeros(numel(kp),numel(wperp));
    Kn = @(x__,y__) kernel(x__,y__);
    % dnjs_ = @(n,j,s) dnjs(n,j,s);
    dnjs_ = @(n,j,s) dnjs_GYAC(n+1,j+1,s+1);
    for id1 = 1:numel(js)
    for id2 = 1:numel(wperp)
        j_ = js(id1);
        v_ = wperp(id2);
        % J0 x Lj
        J0Lj_ = besselj(0,kp*v_)*polyval(LaguerrePoly(j_),v_^2);
        % Taylor exp of J0
        T1_ = (1 - kp.^2*v_^2/4)*polyval(LaguerrePoly(j_),v_^2);
        % Sum_n^Nmax Kn x Ln x Lj
        KLL = zeros(size(kp));
        for n_ = 0:Jmax
            KLL = KLL + Kn(n_,kp).*polyval(LaguerrePoly(n_),v_^2);
        end
        KLL = KLL.*polyval(LaguerrePoly(j_),v_^2);
        % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
        KdL1_ = zeros(size(kp));
        for n_ = 0:Jmax-j_
            tmp_ = 0;
            for s_ = 0:n_+j_
                tmp_ = tmp_ + dnjs_(n_,j_,s_)*polyval(LaguerrePoly(s_),v_^2);
            end
            KdL1_ = KdL1_ + Kn(n_,kp).*tmp_;
        end
        % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
        KdL2_ = zeros(size(kp));
        for n_ = 0:Jmax
            tmp_ = 0;
            for s_ = 0:min(Jmax,n_+j_)
                tmp_ = tmp_ + dnjs_(n_,j_,s_)*polyval(LaguerrePoly(s_),v_^2);
            end
            KdL2_ = KdL2_ + Kn(n_,kp).*tmp_;
        end
        f_  = @(x) max(x);
        kp2 = kp.^2+0.001; intJ0Lj = f_(abs(J0Lj_));
        err_KLL_(id2)  = f_(abs(KLL -J0Lj_)./kp2)./f_(intJ0Lj./kp2);
        err_KdL1_(id2) = f_(abs(KdL1_-J0Lj_)./kp2)./f_(intJ0Lj./kp2);
        err_KdL2_(id2) = f_(abs(KdL2_-J0Lj_)./kp2)./f_(intJ0Lj./kp2);
        err_T1_(id2)   = f_(abs(T1_  -J0Lj_)./kp2)./f_(intJ0Lj./kp2);
        J0Lj(:,id2)    = J0Lj_;
        KdL1(:,id2)    = KdL1_;
        KdL2(:,id2)    = KdL2_;
        T1  (:,id2)    = T1_;
    end
        err_KLL(id1)  = sum(err_KLL_)./numel(wperp);
        err_KdL1(id1) = sum(err_KdL1_)./numel(wperp);
        err_KdL2(id1) = sum(err_KdL2_)./numel(wperp);
        err_T1(id1) = sum(err_T1_)./numel(wperp);
    end
    clr_ = lines(10);
    figure
    plot(js,err_KLL,'-',...
        'DisplayName','$\sum_{n=0}^{J}\mathcal K_n L_n L_j$'); hold on
    plot(js,err_KdL1,':x',...
        'DisplayName','$\sum_{n=0}^{J-j}\mathcal K_n \sum_{s=0}^{n+j} d_{njs} L_s$',...
        'Color',clr_(5,:)); hold on
    plot(js,err_KdL2,':o',...
        'DisplayName','$\sum_{n=0}^{J}\mathcal K_n \sum_{s=0}^{\min(J,n+j)} d_{njs} L_s$',...
        'Color',clr_(2,:)); hold on
    plot(js,err_T1,':s',...
        'DisplayName','$(1 - l_\perp) L_j$',...
        'Color',clr_(4,:)); hold on
    set(gca,'YScale','log');
    xlabel('$j$'); ylabel('$\varepsilon_J$')
    if 1
    %%
        figure
        nc = 25
        subplot(131)
        contourf(kp,wperp,J0Lj',nc); xlabel('$l_{\perp a}$'); ylabel('$w_{\perp a}$')
        clim_ = clim;
        subplot(132)
        contourf(kp,wperp,KdL1', nc); xlabel('$l_{\perp a}$'); ylabel('$w_{\perp a}$')
        clim(clim_);
        % subplot(132)
        subplot(133)
        contourf(kp,wperp,KdL2', nc); xlabel('$l_{\perp a}$'); ylabel('$w_{\perp a}$')
        % contourf(kp,wperp,T1  ,20); xlabel('$k_\perp$'); ylabel('$w_\perp$')
        clim(clim_);
    end
end

