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
% fname = 'CBC_100_20x1x32x30x14_Lv_3_Lw_12_circ.txt';
% fname = 'CBC_100_20x1x32x32x12_Lv_3_Lw_12.txt';
% fname = 'CBC_KT_4_20x1x32x32x12_Lv_3_Lw_12.txt';
% fname = 'CBC_KT_4_20x1x32x64x24_Lv_6_Lw_24.txt';
% fname = 'CBC_KT_5.3_20x1x32x32x12_Lv_3_Lw_12.txt';
% fname = 'CBC_KT_5.3_32x1x48x40x16_Lv_3_Lw_12.txt';
% fname = 'CBC_ky_0.3_20x1x32x32x12_Lv_3_Lw_12.txt';
% fname = 'CBC_ky_0.3_20x1x32x32x12_Lv_3_Lw_12_nuv_1e-3.txt';
% fname = 'CBC_KT_11_20x1x16x24x10_Lv_3_Lw_12.txt';
% fname = 'CBC_KT_11_20x1x32x30x14_Lv_3_Lw_12.txt';
% fname = 'CBC_ky_0.3_20x1x16x24x10_Lv_3_Lw_12_nuv_1e-3.txt';
%----------Shearless CBC
% fname = 'CBC_salpha_s0_nz_24_nv_48_nw_16_adiabe.txt';
% fname = 'CBC_salpha_s0_nz_24_nv_48_nw_16_kine.txt';
% fname = 'CBC_salpha_s0_nz_24_nv_48_nw_16_kine_beta_1e-4.txt';
% fname = 'CBC_miller_s0_nz_24_nv_48_nw_16_adiabe.txt';
% fname = 'CBC_miller_s0_nz_24_nv_48_nw_16_kine.txt';
%----------Shearless pITG
% fname = 'pITG_salpha_s0_nz_24_nv_48_nw_16_adiabe.txt';
% fname = 'pITG_miller_s0_nz_24_nv_48_nw_16_adiabe.txt';
% fname = 'pITG_salpha_s0_nz_24_nv_48_nw_16_kine.txt';
% fname = 'pITG_miller_s0_nz_24_nv_48_nw_16_kine.txt';
%----------Convergence nvpar shearless pITG
% fname = 'pITG_salpha_s0_nz_24_nv_scan_nw_16_adiabe.txt';
% fname = 'pITG_miller_s0_nz_24_nv_scan_nw_16_adiabe.txt';
% fname = 'pITG_salpha_s0_nz_24_nv_scan_nw_16_kine.txt';
% fname = 'pITG_miller_s0_nz_24_nv_scan_nw_16_kine.txt';
% fname = 'pITG_salpha_s0_nz_24_nv_scan_nw_24_adiabe.txt';
% fname = 'pITG_miller_s0_nz_24_nv_scan_nw_24_adiabe.txt';
%----------Convergence nvpar shearless CBC
% fname = 'CBC_salpha_nz_24_nv_scan_nw_16_adiabe.txt';
% fname = 'CBC_miller_nz_24_nv_scan_nw_16_adiabe.txt';
% fname = 'CBC_salpha_nz_24_nv_scan_nw_16_kine.txt';
% fname = 'CBC_miller_nz_24_nv_scan_nw_16_kine.txt';
%---------- CBC
% fname = 'CBC_salpha_nx_8_nz_18_nv_12_nw_8_adiabe.txt';
% fname = 'CBC_salpha_nx_8_nz_18_nv_18_nw_8_adiabe.txt';

% fname = 'CBC_salpha_nx_8_nz_24_nv_8_nw_4_adiabe.txt';
% fname = 'CBC_salpha_nx_8_nz_24_nv_16_nw_8_adiabe.txt';
% fname = 'CBC_salpha_nx_8_nz_24_nv_36_nw_16_adiabe.txt';
fname = 'CBC_salpha_nx_8_nz_24_nv_48_nw_16_adiabe.txt';

% fname = 'kT_5.3_salpha_nx_8_nz_24_nv_60_nw_30_adiabe.txt';

% fname = 'CBC_salpha_nx_8_nz_24_nv_36_nw_16_kine.txt';
% fname = 'CBC_miller_nx_20_nz_32_nv_32_nw_12_adiabe.txt';
% fname = 'CBC_miller_nx_8_nz_24_nv_36_nw_16_adiabe.txt';
% fname = 'CBC_miller_nx_20_nz_32_nv_32_nw_12_kine.txt';
% fname = 'CBC_miller_nx_8_nz_24_nv_36_nw_16_kine.txt';
%----------Convergence nv CBC
% fname = 'CBC_ky_0.3_nv_scan_8x1x24_nw_8_Lv_3_Lw_6.txt';
% fname = 'CBC_ky_0.3_nv_scan_8x1x24_nw_16_Lv_3_Lw_6.txt';
% fname = 'CBC_ky_0.3_nv_scan_8x1x24_nw_24_Lv_3_Lw_6.txt';

data_ = load([path,fname]);

figure
plot(data_(:,2),data_(:,3),'-dk','DisplayName',fname); hold on;
plot(data_(:,2),data_(:,4),'--*k','DisplayName',fname);