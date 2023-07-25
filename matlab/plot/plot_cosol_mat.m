% MAT = results.iCa;
% figure
% suptitle('FCGK,P=6,J=3');
% subplot(221)
%     imagesc(log(abs(MAT))); 
%     title('log abs')
% subplot(222)
%     imagesc(imag(MAT));  colormap(bluewhitered)
%     title('imag'); colorbar
% subplot(223)
%     imagesc(imag(MAT)>0);
%     title('imag$>$0');
% subplot(224)
%     imagesc(imag(MAT)<0);
%     title('imag$<$0');
%     
%     
    %% SGGK
P_ = 20; J_ = 10;
mat_file_name = '/home/ahoffman/HeLaZ/iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5';
SGGK_self = h5read(mat_file_name,'/00000/Caapj/Ciipj');
SGGK_ei = h5read(mat_file_name,'/00000/Ceipj/CeipjF')+h5read(mat_file_name,'/00000/Ceipj/CeipjT');
SGGK_ie = h5read(mat_file_name,'/00000/Ciepj/CiepjF')+h5read(mat_file_name,'/00000/Ciepj/CiepjT');

figure
MAT = 1i*SGGK_self; suptitle('SGGK ii,P=20,J=10, k=0');
% MAT = 1i*SGGK_ei; suptitle('SGGK ei,P=20,J=10');
% MAT = 1i*SGGK_ie; suptitle('SGGK ie,P=20,J=10');
subplot(221)
    imagesc(abs(MAT));
    title('log abs')
subplot(222)
    imagesc(imag(MAT)); colormap(bluewhitered)
    title('imag'); colorbar
subplot(223)
    imagesc(imag(MAT)>0);
    title('imag$>$0');
subplot(224)
    imagesc(imag(MAT)<0);
    title('imag$<$0');

        %% PAGK
P_ = 20; J_ = 10;
mat_file_name = '/home/ahoffman/HeLaZ/iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5';
PAGK_self = h5read(mat_file_name,'/00000/Caapj/Ceepj');
% PAGK_ei = h5read(mat_file_name,'/00000/Ceipj/CeipjF')+h5read(mat_file_name,'/00000/Ceipj/CeipjT');
% PAGK_ie = h5read(mat_file_name,'/00000/Ciepj/CiepjF')+h5read(mat_file_name,'/00000/Ciepj/CiepjT');

figure
MAT = 1i*PAGK_self; suptitle('PAGK ii,P=20,J=10, k=0');

subplot(221)
    imagesc(abs(MAT));
    title('log abs')
subplot(222)
    imagesc(imag(MAT)); colormap(bluewhitered)
    title('imag'); colorbar
subplot(223)
    imagesc(imag(MAT)>0);
    title('imag$>$0');
subplot(224)
    imagesc(imag(MAT)<0);
    title('imag$<$0');

    
    %% FCGK
P_ = 4; J_ = 2;
mat_file_name = '/home/ahoffman/gyacomo/iCa/gk_coulomb_NFLR_12_P_4_J_2_N_50_kpm_4.0.h5';
% mat_file_name = '/home/ahoffman/gyacomo/iCa/LDGK_P6_J3_dk_5e-2_km_2.5_NFLR_20.h5';
% mat_file_name = '/home/ahoffman/gyacomo/iCa/LDGK_P10_J5_dk_5e-2_km_5_NFLR_12.h5';

kp = 0.0;
kp_a =  h5read(mat_file_name,'/coordkperp');
[~,matidx] = min(abs(kp_a-kp));
matidx = sprintf('%5.5i',matidx);
FCGK_self = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
FCGK_ei = h5read(mat_file_name,['/',matidx,'/Ceipj/CeipjF'])+h5read(mat_file_name,'/00000/Ceipj/CeipjT');
FCGK_ie = h5read(mat_file_name,['/',matidx,'/Ciepj/CiepjF'])+h5read(mat_file_name,'/00000/Ciepj/CiepjT');
figure
MAT = 1i*FCGK_self; suptitle(['FCGK ii,P=6,J=3, k=',num2str(kp),' (',matidx,')']);
% MAT = 1i*FCGK_ei; suptitle('FCGK ei,P=20,J=10');
% MAT = 1i*FCGK_ie; suptitle('FCGK ie,P=20,J=10');
subplot(221)
    imagesc(abs(MAT));
    title('log abs')
subplot(222)
    imagesc(imag(MAT)); colormap(bluewhitered)
    title('imag'); colorbar
subplot(223)
    imagesc(imag(MAT)>0);
    title('imag$>$0');
subplot(224)
    imagesc(imag(MAT)<0);
    title('imag$<$0');

%% Eigenvalue spectrum analysis    
if 0
%%
mfns = {...
        '/home/ahoffman/gyacomo/iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5';...
        '/home/ahoffman/gyacomo/iCa/LDGK_P6_J3_dk_5e-2_km_2.5_NFLR_20.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_coulomb_NFLR_12_P_4_J_2_N_50_kpm_4.0.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_landau_P10_J5_dk_5e-2_km_2.0_NFLR_12.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_landau_P11_J7_dk_5e-2_km_2.0_NFLR_16.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_landauii_P16_J9_dk_5e-2_km_2.0_NFLR_8.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_landau_P16_J9_dk_5e-2_km_2.0_NFLR_8.h5';...
%         '/home/ahoffman/gyacomo/iCa/gk.hacked_sugama_P_10_J_5_N_150_kpm_8.0.h5';...
%         '/home/ahoffman/gyacomo/iCa/gk.hacked_sugama_P_4_J_2_N_75_kpm_5.0.h5';...
        };
CONAME_A = {...
    'SG 20 10';...
    'PA 20 10';...
    'LD 6  3 NFLR 20'; ...
    'FC 4 2 NFLR 12'; ...
    'LD 10 5 NFLR 12'; ...
    'LD 11 7 NFLR 16'; ...
    'LDii 16 9 NFLR 8'; ...
    'LD   16 9 NFLR 8'; ...
%     'Hacked SG A';...
%     'Hacked SG B';...
    };
figure
for j_ = 1:numel(mfns)
    mat_file_name = mfns{j_}; disp(mat_file_name);
    kp_a =  h5read(mat_file_name,'/coordkperp');
%     kp_a = kp_a(kp_a<=3);
    gmax = zeros(size(kp_a));
    wmax = zeros(size(kp_a));
   for idx_ = 0:numel(kp_a)-1
        matidx = sprintf('%5.5i',idx_);
        MAT = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
        gmax(idx_+1) = max((real(eig(MAT)))); 
        wmax(idx_+1) = max(abs(imag(eig(MAT)))); 
   end
   subplot(121)
    plot(kp_a,gmax,'DisplayName',CONAME_A{j_}); hold on;
   subplot(122)
    plot(kp_a,wmax,'DisplayName',CONAME_A{j_}); hold on;
end
   subplot(121)
legend('show'); grid on;
ylim([0,100]);
xlabel('$k_\perp$'); ylabel('$\gamma_{max}$ from Eig(iCa)')
   subplot(122)
legend('show'); grid on;
xlabel('$k_\perp$'); ylabel('$\omega_{max}$ from Eig(iCa)')
end

%% Van Kampen plot
if 0
%%
kperp= 1.5;
mfns = {...
        '/home/ahoffman/gyacomo/iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_coulomb_NFLR_6_P_4_J_2_N_50_kpm_4.0.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_coulomb_NFLR_12_P_4_J_2_N_50_kpm_4.0.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_landau_P10_J5_dk_5e-2_km_2.0_NFLR_12.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_landauii_P16_J9_dk_5e-2_km_2.0_NFLR_8.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_sugama_tau1e-3_P11_J7_dk_5e-2_km_5.0_NFLR_12.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_sugama_tau1e-3_P4_J2_dk_5e-2_km_5.0_NFLR_5.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_sugama_tau1e-2_P4_J2_dk_5e-2_km_5.0_NFLR_5.h5';...
        };
CONAME_A = {'SG 20 10';...
    'PA 20 10';...
    'FC 4 2 NFLR 6';...
    'FC 4 2 NFLR 12';...
    'LD 10 5 NFLR 12';...
    'LD 16 9 NFLR 16';...
    'SG 11 7 NFLR 12, tau 1e-3'; ...
    'SG 4 2 NFLR 5, tau 1e-3'; ...
    'SG 4 2 NFLR 5, tau 1e-2'; ...
    };
TAU_A = [1;...
    1;...
    1;...
    1;...
    1;...
    1;...
    1e-3;...
    1e-3;...
    1e-2;...
    ];
grow = {};
puls = {};
for j_ = 1:numel(mfns)
    mat_file_name = mfns{j_}; disp(mat_file_name);
    kp_a          =  h5read(mat_file_name,'/coordkperp');
    [~,idx_]       = min(abs(kp_a-kperp));
    matidx        = sprintf('%5.5i',idx_);
    MAT           = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
    grow{j_}  = real(eig(MAT))*TAU_A(j_); 
    puls{j_}  = imag(eig(MAT))*TAU_A(j_); 
end

figure
for j_ = 1:numel(mfns)
%    plot(puls{j_}, grow{j_},'o','DisplayName',CONAME_A{j_}); hold on;
   plot(grow{j_},'o','DisplayName',CONAME_A{j_}); hold on;
end
legend('show'); grid on; title(['$k_\perp=$',num2str(kperp)]);
xlabel('$\omega$ from Eig(iCa)'); ylabel('$\gamma$ from Eig(iCa)')
end