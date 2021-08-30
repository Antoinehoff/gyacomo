if 0
%% Kernels
kmax=7;
nmax=6;
kx = linspace(0,kmax,100);


figure
for n_ = 0:nmax
    plot(kx,kernel(n_,kx),'DisplayName',['$\mathcal{K}_{',num2str(n_),'}$']);hold on;
end
ylim_ = ylim;
plot(kx(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
plot(kx,J0,'-r','DisplayName','$J_0$');
legend('show')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
%% Bessels and approx
vperp = linspace(0,1.5,4);
nmax1=5;
nmax2=10;
kmax=7;
figure
for i = 1:4
subplot(2,2,i)
    v_ = vperp(i);
    kx = linspace(0,kmax,100);

    J0 = besselj(0,kx*v_);
    A1 = 1 - kx.^2*v_^2/4;
    K1 = zeros(size(kx));
    K2 = zeros(size(kx));
    for n_ = 0:nmax1
        K1 = K1 + kernel(n_,kx).*polyval(LaguerrePoly(n_),v_^2);
    end
    for n_ = 0:nmax2
        K2 = K2 + kernel(n_,kx).*polyval(LaguerrePoly(n_),v_^2);
    end
    plot(kx,J0,'-k','DisplayName','$J_0(kv)$'); hold on;
    plot(kx,A1,'-r','DisplayName','$1 - k^2 v^2/4$');
    plot(kx,K1,'--b','DisplayName',['$\sum_{n=0}^{',num2str(nmax1),'}\mathcal K_n(k) L^n(v)$']);
    plot(kx,K2,'-b','DisplayName',['$\sum_{n=0}^{',num2str(nmax2),'}\mathcal K_n(k) L^n(v)$']);
    ylim_ = [-0.5, 1.0];
    plot(kx(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
    ylim(ylim_); xlabel('$k$')
    legend('show'); grid on; title(['$v = ',num2str(v_),'$'])
end
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%% from Sum_n Kn x Ln x Lj to Sum_n Kn Sum_s dnjs x Ls
vperp = [5];
kx    = linspace(0,10,200);
Jmax  = 15;
j     = 10;
fig = figure; set(gcf, 'Position',  [100, 100, numel(vperp)*300, 250])
%     suptitle(['$J_{max}=',num2str(Jmax),', j=',num2str(j),'$'])
Kn = @(x__,y__) kernel(x__,y__);
for i = 1:numel(vperp)
subplot(1,numel(vperp),i)
    v_ = vperp(i);
    % J0 x Lj
    J0Lj = besselj(0,kx*v_)*polyval(LaguerrePoly(j),v_^2);
    % Taylor exp of J0
    A1 = (1 - kx.^2*v_^2/4)*polyval(LaguerrePoly(j),v_^2);
    % Sum_n^Nmax Kn x Ln x Lj
    KLL = zeros(size(kx));
    for n_ = 0:Jmax
        KLL = KLL + Kn(n_,kx).*polyval(LaguerrePoly(n_),v_^2);
    end
    KLL = KLL.*polyval(LaguerrePoly(j),v_^2);
    % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
    KdL1 = zeros(size(kx));
    for n_ = 0:Jmax-j
        tmp_ = 0;
        for s_ = 0:n_+j
            tmp_ = tmp_ + dnjs(n_,j,s_)*polyval(LaguerrePoly(s_),v_^2);
        end
        KdL1 = KdL1 + Kn(n_,kx).*tmp_;
    end
    % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
    KdL2 = zeros(size(kx));
    for n_ = 0:Jmax
        tmp_ = 0;
        for s_ = 0:min(Jmax,n_+j)
            tmp_ = tmp_ + dnjs(n_,j,s_)*polyval(LaguerrePoly(s_),v_^2);
        end
        KdL2 = KdL2 + Kn(n_,kx).*tmp_;
    end
    
    % plots
    plot(kx,J0Lj,'-k','DisplayName','$J_0 L_j$'); hold on;
    plot(kx,KLL,'-','DisplayName',['$\sum_{n=0}^{J}\mathcal K_n L^n L_j$']);
    plot(kx,KdL1,'-.','DisplayName',['$\sum_{n=0}^{J-j}\mathcal K_n \sum_{s=0}^{n+j} d_{njs} L^s$']);
    ylim_ = ylim;
    plot(kx,A1,'-','DisplayName','$(1 - k^2 v^2/4) L_j$'); hold on;
    plot(kx,KdL2,'--','DisplayName',['$\sum_{n=0}^{J}\mathcal K_n \sum_{s=0}^{\min(J,n+j)} d_{njs} L^s$']);
%     plot(kx(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
    ylim(ylim_); 
    xlabel('$k_{\perp}$')
    legend('show'); 
    grid on; title(['$v = ',num2str(v_),'$',' $J_{max}=',num2str(Jmax),', j=',num2str(j),'$'])
end
FIGNAME = ['/home/ahoffman/Dropbox/Applications/Overleaf/Hermite-Laguerre Z-pinch (HeLaZ Doc.)/figures/J0Lj_approx_J','_',num2str(Jmax),'_j_',num2str(j),'.eps'];
% saveas(fig,FIGNAME,'epsc');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
%% Convergence analysis
kx    = linspace(0,10,200);
v_A   = linspace(0,1,10); 
J_A   = 1:15;
[YY,XX] = meshgrid(v_A,J_A);
klimLL = zeros(size(XX));
klimdL1= zeros(size(XX));
Kn = @(x__,y__) kernel(x__,y__);
fig = figure; set(gcf, 'Position',  [100, 100, 500, 400]);
for j = 5
for ij_ = 1:numel(J_A)
    iv_ = 1;
    for v_ = v_A
    eps   = abs(0.01*polyval(LaguerrePoly(j),v_^2));
    Jmax = J_A(ij_);
    % J0 x Lj
    J0Lj = besselj(0,kx*v_)*polyval(LaguerrePoly(j),v_^2);
    % Taylor exp of J0
    A1 = (1 - kx.^2*v_^2/4)*polyval(LaguerrePoly(j),v_^2);
    % Sum_n^Nmax Kn x Ln x Lj
    KLL = zeros(size(kx));
    for n_ = 0:Jmax
        KLL = KLL + Kn(n_,kx).*polyval(LaguerrePoly(n_),v_^2);
    end
    KLL = KLL.*polyval(LaguerrePoly(j),v_^2);
    % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
    KdL1 = zeros(size(kx));
    for n_ = 0:Jmax-j
        tmp_ = 0;
        for s_ = 0:n_+j
            tmp_ = tmp_ + dnjs(n_,j,s_)*polyval(LaguerrePoly(s_),v_^2);
        end
        KdL1 = KdL1 + Kn(n_,kx).*tmp_;
    end
    % Sum_n^Nmax Kn Sum_s^Smax dnjs x Ls
    KdL2 = zeros(size(kx));
    for n_ = 0:Jmax
        tmp_ = 0;
        for s_ = 0:min(Jmax,n_+j)
            tmp_ = tmp_ + dnjs(n_,j,s_)*polyval(LaguerrePoly(s_),v_^2);
        end
        KdL2 = KdL2 + Kn(n_,kx).*tmp_;
    end
    % errors
    err_     = abs(J0Lj-KLL);
    idx_ = 1;
    while(err_(idx_)< eps && idx_ < numel(err_))
        idx_ = idx_ + 1;
    end
    klimLL(ij_,iv_) = kx(idx_);
    err_     = abs(J0Lj-KdL1);
    idx_ = 1;
    while(err_(idx_)< eps && idx_ < numel(err_))
        idx_ = idx_ + 1;
    end
    klimdL1(ij_,iv_) = kx(idx_);
    iv_ = iv_ + 1;
    end
end
% plot
% plot(JmaxA,klimLL,'o-'); hold on
% plot(J_A,klimdL1,'o-'); hold on
end
surf(XX,YY,klimdL1)
xlabel('$J_{max}$'); ylabel('$v_{\perp}$');
title(['$k$ s.t. $\epsilon=1\%$, $j=',num2str(j),'$'])
ylim([0,1]); grid on;
end
if 0
%% Test Lj Ln = sum_s dnjs Ls
j = 10;
n = 10;
smax = n+j;
vperp = linspace(0,5,20);
LjLn        = zeros(size(vperp));
dnjsLs      = zeros(size(vperp));
dnjsLs_bin  = zeros(size(vperp));
dnjsLs_stir = zeros(size(vperp));

i=1;
for x_ = vperp
    LjLn(i) = polyval(LaguerrePoly(j),x_)*polyval(LaguerrePoly(n),x_);
    dnjsLs(i)      = 0;
    dnjsLs_bin(i)  = 0;
    dnjsLs_stir(i) = 0;
    for s_ = 0:smax
        dnjsLs(i)      = dnjsLs(i)      + dnjs(n,j,s_)*polyval(LaguerrePoly(s_),x_);
        dnjsLs_bin(i)  = dnjsLs_bin(i)  + dnjs_fact(n,j,s_)*polyval(LaguerrePoly(s_),x_);
        dnjsLs_stir(i) = dnjsLs_stir(i) + dnjs_stir(n,j,s_)*polyval(LaguerrePoly(s_),x_);
    end
    i = i+1;
end

figure
plot(vperp,LjLn); hold on;
plot(vperp,dnjsLs,'d')
plot(vperp,dnjsLs_bin,'o')
plot(vperp,dnjsLs_stir/dnjsLs_stir(1),'x')
% We see that starting arround j = 18, n = 0, the stirling formula is the only
% approximation that does not diverge from the target function. However it
% performs badly for non 0 n...
end



