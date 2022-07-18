%% gyroLES plot
% genepath = '/misc/gene_results';
% % simpath  = '/NL_Zpinch_Kn_1.7_eta_0.25_nuSG_1e-1_gyroLES/';
% simpath  = '/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-2_gyroLES/';
% filename = 'GyroLES_ions.dat';
% toload = [genepath simpath filename];
% data = readtable(toload);
% t = data.Var1;
% mu_x = data.Var2;
% mu_y = data.Var3;
% figure;
% plot(t,mu_x); hold on
% plot(t,mu_y);
% legend('mu_x','mu_y'); xlabel('time'); ylabel('Hyper-diff');
% title('gyroLES results for H&P fig. 2c')

%% Plot linear results
% path = '/home/ahoffman/gene/linear_zpinch_results/';
% fname ='tmp.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.45_nuSG_0.047_32x16.txt';
% fname ='GENE_LIN_Kn_1.5_KT_0.325_nuSG_0.047_32x16.txt';
% fname ='GENE_LIN_Kn_1.5_KT_0.325_nuSG_0.0235_32x16.txt';
% fname ='GENE_LIN_Kn_1.5_KT_0.325_nuSG_0.0047_32x16.txt';
% fname ='GENE_LIN_Kn_1.5_KT_0.325_nuSG_0.0047_64x32.txt';
% fname ='GENE_LIN_Kn_1.4_KT_0.35_nuSG_0.0047_64x32.txt';
% fname ='GENE_LIN_Kn_1.6_KT_0.4_nuSG_0.0047_64x32.txt';
% fname ='GENE_LIN_Kn_1.6_KT_0.4_nuSGDK_0.0047_64x32.txt';
% fname ='GENE_LIN_Kn_1.6_KT_0.4_nuLDDK_0.0047_64x32.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.4_nuLDDK_0.0047_64x32.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.45_nu_0_32x16.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.45_nuSG_0.0235_64x32.txt';
% fname ='GENE_LIN_Kn_1.9_KT_0.475_nuSG_0.047_64x32.txt';
% fname ='GENE_LIN_Kn_1.9_KT_0.475_nuSG_0.235_64x32.txt';
% fname ='GENE_LIN_Kn_1.7_KT_0.425_nuSG_0.235_64x32.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.45_nuSGDK_0.047_32x16.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.475_nu_0_mu_5e-2.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.45_nuSGDK_0.0235_64x32.txt';
% fname ='GENE_LIN_Kn_1.8_KT_0.45_nuSGDK_0.0235_32x16.txt';
% fname ='GENE_LIN_Kn_2.0_KT_0.5_nu_0_32x16.txt';
% fname ='GENE_LIN_Kn_2.0_KT_0.5_nuSGDK_0.0235_32x16.txt';
% fname ='GENE_LIN_Kn_1.6_KT_0.4_nu_0_32x16.txt';
% fname ='GENE_LIN_Kn_2.5_KT_0.625_nu_0_32x16.txt';
path = '/home/ahoffman/gene/linear_CBC_results/';
fname = 'CBC_linear.txt';
data_ = load([path,fname]);

figure
plot(data_(:,2),data_(:,3),'-dk','DisplayName','GENE');