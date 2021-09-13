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
P_ = 6; J_ = 3;
mat_file_name = '/home/ahoffman/HeLaZ/iCa/gk_coulomb_P_6_J_3_N_150_kpm_8.0_NFLR_8.h5';
kp = 1.0; matidx = round(kp/(8/150));
% matidx = sprintf('%5.5i',matidx);
% mat_file_name = '/home/ahoffman/HeLaZ/iCa/gk_coulomb_P_6_J_3_N_1_kpm_4.0_NFLR_20_m_0.h5';
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
    
%% Single eigenvalue analysis

% mat_file_name = '/home/ahoffman/cosolver/gk.coulomb_6_3_NFLR_8_kp_4_bjf/scanfiles_00000/self.0.h5';
mat_file_name = '/home/ahoffman/cosolver/gk.coulomb_6_3_NFLR_8_kp_4_new_Tljpmf_fort_bjf/scanfiles_00000/self.0.h5';
matidx = 74;

matidx = sprintf('%5.5i',matidx);disp(matidx);

% MAT = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
MAT = h5read(mat_file_name,['/Caapj/Ciipj']);

gmax = max(real(eig(MAT)));

wmax = max(imag(eig(MAT)));
figure
imagesc((MAT)); colormap(bluewhitered)
title(['$\gamma_{max}=',num2str(gmax),'$'])
%% Eigenvalue spectrum analysis    
if 0
%%
mfns = {'/home/ahoffman/HeLaZ/iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5',...
        '/home/ahoffman/HeLaZ/iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5',...
        '/home/ahoffman/HeLaZ/iCa/gk_coulomb_P_10_J_5_N_20_kpm_8.0_NFLR_0.h5',...
        '/home/ahoffman/HeLaZ/iCa/gk_coulomb_P_2_J_1_N_20_kpm_8.0_NFLR_5.h5',...
        '/home/ahoffman/HeLaZ/iCa/gk_coulomb_P_6_J_3_N_150_kpm_8.0_NFLR_0.h5',...
        '/home/ahoffman/HeLaZ/iCa/gk_coulomb_P_6_J_3_N_150_kpm_8.0_NFLR_8.h5',...
        };
CONAME_A = {'SG 20 10','PA 20 10', 'FC 10 5 NFLR 0', 'FC 2 1 NFLR 5', 'FC 6 3 NFLR 0', 'FC 6 3 NFLR 8'};
figure
for j_ = 1:numel(mfns)
    mat_file_name = mfns{j_}; disp(mat_file_name);
    kp_a =  h5read(mat_file_name,'/coordkperp');
    gmax = zeros(size(kp_a));
    wmax = zeros(size(kp_a));
   for idx_ = 0:numel(kp_a)-1
        matidx = sprintf('%5.5i',idx_);
        MAT = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
        gmax(idx_+1) = max(real(eig(MAT))); 
        wmax(idx_+1) = max(imag(eig(MAT))); 
   end
%    subplot(121)
    plot(kp_a,gmax,'DisplayName',CONAME_A{j_}); hold on;
%    subplot(122)
%     plot(kp_a,wmax,'DisplayName',CONAME_A{j_}); hold on;
end
%    subplot(121)
legend('show'); grid on;
xlabel('$k_\perp$'); ylabel('$\gamma_{max}$ from Eig(iCa)')
%    subplot(122)
% legend('show'); grid on;
% xlabel('$k_\perp$'); ylabel('$\omega_{max}$ from Eig(iCa)')
end
