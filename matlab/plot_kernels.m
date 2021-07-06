%% Kernels
kmax=7;
nmax=6;
kr_ = linspace(0,kmax,100);


figure
for n_ = 0:nmax
    plot(kr_,kernel(n_,kr_),'DisplayName',['$\mathcal{K}_{',num2str(n_),'}$']);hold on;
end
ylim_ = ylim;
plot(kr_(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
plot(kr_,J0,'-r','DisplayName','$J_0$');
legend('show')

%% Bessels and approx
vperp = linspace(0,1.5,4);
nmax1=5;
nmax2=10;
kmax=7;
figure
for i = 1:4
subplot(2,2,i)
    v_ = vperp(i);
    kr_ = linspace(0,kmax,100);

    J0 = besselj(0,kr_*v_);
    A1 = 1 - kr_.^2*v_^2/4;
    K1 = zeros(size(kr_));
    K2 = zeros(size(kr_));
    for n_ = 0:nmax1
        K1 = K1 + kernel(n_,kr_).*polyval(LaguerrePoly(n_),v_^2);
    end
    for n_ = 0:nmax2
        K2 = K2 + kernel(n_,kr_).*polyval(LaguerrePoly(n_),v_^2);
    end
    plot(kr_,J0,'-k','DisplayName','$J_0(kv)$'); hold on;
    plot(kr_,A1,'-r','DisplayName','$1 - k^2 v^2/4$');
    plot(kr_,K1,'--b','DisplayName',['$\sum_{n=0}^{',num2str(nmax1),'}\mathcal K_n(k) L^n(v)$']);
    plot(kr_,K2,'-b','DisplayName',['$\sum_{n=0}^{',num2str(nmax2),'}\mathcal K_n(k) L^n(v)$']);
    ylim_ = [-0.5, 1.0];
    plot(kr_(end)*[2/3 2/3],ylim_,'--k','DisplayName','AA');
    ylim(ylim_); xlabel('$k$')
    legend('show'); grid on; title(['$v = ',num2str(v_),'$'])
end
%%