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

    
    %% SGGK
P_ = 6; J_ = 3;
mat_file_name = '/home/ahoffman/HeLaZ/iCa/gk_coulomb_P_6_J_3_N_150_kpm_8.0.h5';
FCGK_self = h5read(mat_file_name,'/00000/Caapj/Ciipj');
FCGK_ei = h5read(mat_file_name,'/00000/Ceipj/CeipjF')+h5read(mat_file_name,'/00000/Ceipj/CeipjT');
FCGK_ie = h5read(mat_file_name,'/00000/Ciepj/CiepjF')+h5read(mat_file_name,'/00000/Ciepj/CiepjT');

figure
MAT = 1i*FCGK_self; suptitle('FCGK ii,P=20,J=10, k=0');
% MAT = 1i*SGGK_ei; suptitle('FCGK ei,P=20,J=10');
% MAT = 1i*SGGK_ie; suptitle('FCGK ie,P=20,J=10');
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
    