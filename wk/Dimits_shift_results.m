%% Nu = 1e-2
KN        = [   1.5    1.6    1.7    1.8    1.9    2.0    2.5];
Ginf_LDGK = [5.5e-2 1.3e-1 2.7e-1 4.6e-1 1.0e+0 0.0e+0 0.0e+0];
conv_LDGK = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

Ginf_SGGK = [1.2e-3 2.8e-3 1.7e-2 8.9e-2 3.0e+0 0.0e+0 0.0e+0];
conv_SGGK = [9.0e-4 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];


Ginf_DGGK = [1.1e-4 5.6e-4 9.0e-4 6.5e-3 7.0e+0 0.0e+0 0.0e+0];
conv_DG   = [1.6e-4 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

figure
semilogy(KN,Ginf_LDGK,'d--g','DisplayName','Landau'); hold on;


semilogy(KN,Ginf_SGGK,'o-b','DisplayName','Sugama'); hold on;
semilogy(KN,conv_SGGK,'x-k','DisplayName','conv'); hold on;


semilogy(KN,Ginf_DGGK,'^-r','DisplayName','Dougherty'); hold on;
semilogy(KN,conv_DG,'x-k','DisplayName','conv'); hold on;

title('$\nu=0.01$');
xlabel('Drive');
ylabel('Transport');

%% Nu = 5e-2
KN      = [   1.5    1.6    1.7    1.8    1.9];
Ginf_LDGK = [0.0e+0 0.0e+0 0.0e+0 8.3e-1 0.0e+0];
conv_LDGK = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

Ginf_SGGK = [0.0e+0 0.0e+0 0.0e+0 2.0e-1 0.0e+0];
conv_SGGK = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

Ginf_SGDK = [0.0e+0 0.0e+0 0.0e+0 5.6e-3 0.0e+0];
conv_SGDK = [0.0e+0 0.0e+0 0.0e+0 2.0e-3 0.0e+0];

Ginf_DGGK = [0.0e+0 0.0e+0 0.0e+0 1.5e-2 0.0e+0];
figure
semilogy(KN,Ginf_LDGK,'d--g','DisplayName','Landau GK'); hold on;
semilogy(KN,Ginf_SGGK,'o-b','DisplayName','Sugama GK'); hold on;
semilogy(KN,Ginf_SGDK,'o-c','DisplayName','Sugama DK'); hold on;
semilogy(KN,Ginf_DGGK,'^-r','DisplayName','Dougherty GK'); hold on;
semilogy(KN,conv_LDGK,'x--g','DisplayName','changing params'); hold on;
semilogy(KN,conv_SGGK,'x--b','DisplayName','changing params'); hold on;
semilogy(KN,conv_SGDK,'x--c','DisplayName','changing params'); hold on;
title('$\nu=0.05$');
xlabel('Drive');
ylabel('Transport');


%% Nu = 1e-1
KN      = [   1.5    1.6    1.7    1.8    1.9      2.0    2.5];
Ginf_LDGK = [2.5e-2 2.5e-1 6.3e-1 1.3e+0 2.3e+0 4.3e+0 0.0e+0];
conv_LDGK = [0.0e+0 2.5e-1 0.0e+0 1.2e+0 0.0e+0 0.0e+0 0.0e+0];
Ginf_SGGK = [0.0e-9 1.3e-2 9.1e-2 3.0e-1 7.8e-1 1.4e+0 3.5e+1];
Ginf_DGGK = [2.0e-4 6.4e-3 3.6e-2 3.1e-1 0.0e+0 3.0e-1 3.3e+1];

figure
semilogy(KN,Ginf_LDGK,'d--g','DisplayName','Landau'); hold on;
semilogy(KN,Ginf_SGGK,'o-b','DisplayName','Sugama'); hold on;
semilogy(KN,Ginf_DGGK,'^-r','DisplayName','Dougherty'); hold on;
semilogy(KN,conv_LDGK,'x--k','DisplayName','changing params'); hold on;
title('$\nu=0.1$');
xlabel('Drive');
ylabel('Transport');


%% Nu = 5e-1
KN      = [   1.5    1.6    1.7    1.8    1.9];
Ginf_LDGK = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
conv_LDGK = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
Ginf_SGGK = [0.0e-9 3.9e-2 2.7e-1 6.6e-1 1.6e+0];
conv_SGGK = [0.0e-9 0.0e+0 3.0e-1 0.0e+0 1.6e+0];
Ginf_DGGK = [0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];

figure
semilogy(KN,Ginf_LDGK,'d--g','DisplayName','Landau'); hold on;
semilogy(KN,Ginf_SGGK,'o-b','DisplayName','Sugama'); hold on;
semilogy(KN,Ginf_DGGK,'^-r','DisplayName','Dougherty'); hold on;
semilogy(KN,conv_LDGK,'x--k','DisplayName','changing params'); hold on;
semilogy(KN,conv_SGGK,'x--k','DisplayName','changing params'); hold on;
title('$\nu=0.5$');
xlabel('Drive');
ylabel('Transport');

%% Primary mode stability region (confirmed with GENE)

stable = [1.5 0.1; 1.5 0.05];
unstable = [1.5 0.01];

figure
semilogy(stable(:,1),stable(:,2),'og'); hold on;
semilogy(unstable(:,1),unstable(:,2),'or');
