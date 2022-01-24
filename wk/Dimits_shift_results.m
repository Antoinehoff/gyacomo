%% Nu = 1e-2
KN      = [   1.5    1.6    1.7    1.8    1.9];
Ginf_LD = [5.5e-2 1.3e-1 2.7e-1 4.6e-1 1.0e+0];
conv_LD = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

Ginf_SG = [1.2e-3 2.8e-3 1.7e-2 8.9e-2 3.0e+0];
conv_SG = [9.0e-4 0.0e+0 0.0e+0 0.0e+0 0.0e+0];


Ginf_DG = [1.1e-4 5.6e-4 9.0e-4 6.5e-3 7.0e+0];
conv_DG = [1.6e-4 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

figure
semilogy(KN,Ginf_LD,'d--g','DisplayName','Landau'); hold on;


semilogy(KN,Ginf_SG,'o-b','DisplayName','Sugama'); hold on;
semilogy(KN,conv_SG,'x-k','DisplayName','conv'); hold on;


semilogy(KN,Ginf_DG,'^-r','DisplayName','Dougherty'); hold on;
semilogy(KN,conv_DG,'x-k','DisplayName','conv'); hold on;

title('$\nu=0.01$');
xlabel('Drive');
ylabel('Transport');

%% Nu = 5e-2
KN      = [   1.5    1.6    1.7    1.8    1.9];
Ginf_LD = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
conv_LD = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
Ginf_SG = [0.0e+0 0.0e+0 0.0e+0 2.0e-1 0.0e+0];
conv_SG = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
Ginf_DG = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
figure
semilogy(KN,Ginf_LD,'d--g','DisplayName','Landau'); hold on;
semilogy(KN,Ginf_SG,'o-b','DisplayName','Sugama'); hold on;
semilogy(KN,Ginf_DG,'^-r','DisplayName','Dougherty'); hold on;
semilogy(KN,conv_LD,'x--k','DisplayName','changing params'); hold on;
semilogy(KN,conv_SG,'x--k','DisplayName','changing params'); hold on;
title('$\nu=0.05$');
xlabel('Drive');
ylabel('Transport');


%% Nu = 1e-1
KN      = [   1.5    1.6    1.7    1.8    1.9];
Ginf_LD = [2.5e-2 2.5e-1 6.3e-1 1.3e+0 2.3e+0];
conv_LD = [0.0e+0 2.5e-1 0.0e+0 1.2e+0 0.0e+0];
Ginf_SG = [0.0e-9 1.3e-2 9.1e-2 3.0e-1 7.8e-1];
Ginf_DG = [2.0e-4 6.4e-3 3.6e-2 3.1e-1 0.0e+0];

figure
semilogy(KN,Ginf_LD,'d--g','DisplayName','Landau'); hold on;
semilogy(KN,Ginf_SG,'o-b','DisplayName','Sugama'); hold on;
semilogy(KN,Ginf_DG,'^-r','DisplayName','Dougherty'); hold on;
semilogy(KN,conv_LD,'x--k','DisplayName','changing params'); hold on;
title('$\nu=0.1$');
xlabel('Drive');
ylabel('Transport');


%% Nu = 5e-1
KN      = [   1.5    1.6    1.7    1.8    1.9];
Ginf_LD = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
conv_LD = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
Ginf_SG = [0.0e-9 3.9e-2 2.7e-1 6.6e-1 1.6e+0];
conv_SG = [0.0e-9 0.0e+0 3.0e-1 0.0e+0 1.6e+0];
Ginf_DG = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

figure
semilogy(KN,Ginf_LD,'d--g','DisplayName','Landau'); hold on;
semilogy(KN,Ginf_SG,'o-b','DisplayName','Sugama'); hold on;
semilogy(KN,Ginf_DG,'^-r','DisplayName','Dougherty'); hold on;
semilogy(KN,conv_LD,'x--k','DisplayName','changing params'); hold on;
semilogy(KN,conv_SG,'x--k','DisplayName','changing params'); hold on;
title('$\nu=0.5$');
xlabel('Drive');
ylabel('Transport');

%% Primary mode stability region (confirmed with GENE)

stable = [1.5 0.1; 1.5 0.05];
unstable = [1.5 0.01];

figure
semilogy(stable(:,1),stable(:,2),'og'); hold on;
semilogy(unstable(:,1),unstable(:,2),'or');
