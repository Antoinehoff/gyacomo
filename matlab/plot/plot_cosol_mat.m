%% Coeff trajectory
% matfilename = 'gk_landau_P10_J5_dk_5e-2_km_2.0_NFLR_12.h5'; JMAX = 5;
% matfilename = 'gk_sugama_P_20_J_10_N_150_kpm_8.0.h5'; JMAX = 10;
% matfilename = 'gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5'; JMAX = 10;
matfilename = 'gk_landau_P16_J9_dk_5e-2_km_2.0_NFLR_8.h5'; JMAX = 10;
P_ = 10; J_ = P_/2;
n = 1:((P_+1)*(J_+1));
for p_ = 0:P_
    n((0:J_)+p_*(J_+1)+1) = (0:J_)+p_*(JMAX+1);
end
mat_file_name = ['/home/ahoffman/gyacomo/iCa/',matfilename];
kp_a   = h5read(mat_file_name,'/coordkperp');
coeff  = kp_a*0;
gmax = kp_a*0;

p1 = 0; j1 = 0;
p2 = 0; j2 = 0;
for ik = 1:numel(kp_a)
    matidx = sprintf('%5.5i',ik-1);
    M_self = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
    MAT  = M_self(n+1,n+1); 
    ipj1 = n(j1+p1*(J_+1)+1);
    ipj2 = n(j2+p2*(J_+1)+1);
    coeff(ik)  = MAT(ipj1+1,ipj2+1);
    gmax(ik) = -max((real(eig(MAT))));
    
end
figure
plot(kp_a,gmax);
    %% SGGK
P_ = 10; J_ = 5;
n = 1:((P_+1)*(J_+1));
for p_ = 0:P_
    n((0:J_)+p_*(J_+1)+1) = (0:J_)+p_*(10+1);
end
kp = 5.0;
kp_a =  h5read(mat_file_name,'/coordkperp');
[~,matidx] = min(abs(kp_a-kp));
matidx = sprintf('%5.5i',matidx);
mat_file_name = '/home/ahoffman/gyacomo/iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5';
M_self = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
figure
MAT = M_self; %suptitle('SGGK ii,P=20,J=10, k=0');
imagesc((MAT(n+1,n+1)));
title('Sugama')
xlim('tight')
ylim('tight')
axis equal
% clim(1*[-1 1])
clim auto
colormap(bluewhitered)

        %% PAGK
P_ = 10; J_ = 5;
n = 1:((P_+1)*(J_+1));
for p_ = 0:P_
    n((0:J_)+p_*(J_+1)+1) = (0:J_)+p_*(10+1);
end
kp = 1.0;
kp_a =  h5read(mat_file_name,'/coordkperp');
[~,matidx] = min(abs(kp_a-kp));
matidx = sprintf('%5.5i',matidx);
mat_file_name = '/home/ahoffman/gyacomo/iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5';
PAGK_self = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
figure
PA_MAT = PAGK_self;
    imagesc((PA_MAT(n+1,n+1)));
    title('Lorentz')
    xlim('tight')
    ylim('tight')
    axis equal
    % clim(1*[-1 1])
    clim auto
    colormap(bluewhitered)

    
    %% FCGK
P_ = 10; J_ = 5;
n = 1:((P_+1)*(J_+1));
for p_ = 0:P_
    n((0:J_)+p_*(J_+1)+1) = (0:J_)+p_*(5+1);
end
mat_file_name = '/home/ahoffman/gyacomo/iCa/gk_landau_P10_J5_dk_5e-2_km_2.0_NFLR_12.h5';
kp = 1.0;
kp_a =  h5read(mat_file_name,'/coordkperp');
[~,matidx] = min(abs(kp_a-kp));
matidx = sprintf('%5.5i',matidx);
FCGK_self = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
figure
FC_MAT = FCGK_self;
% MAT = 1i*FCGK_ei; suptitle('FCGK ei,P=20,J=10');
% MAT = 1i*FCGK_ie; suptitle('FCGK ie,P=20,J=10');
    imagesc((FC_MAT(n+1,n+1)));
    title('Landau')
    xlim('tight')
    ylim('tight')
    axis equal
    % clim(1*[-1 1])
    clim auto
   colormap(bluewhitered)

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
        '/home/ahoffman/gyacomo/iCa/gk_landau_P12_J6_dk_5e-2_km_1.0_NFLR_16.h5';...
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
    'LD 12 6 NFLR 16'; ...
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
    c44  = zeros(size(kp_a));
   for idx_ = 0:numel(kp_a)-1
        matidx = sprintf('%5.5i',idx_);
        FC_MAT = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
        gmax(idx_+1) = max((real(eig(FC_MAT)))); 
        wmax(idx_+1) = max((imag(eig(FC_MAT)))); 
        c44 (idx_+1) = FC_MAT(4,4); 
   end
   subplot(121)
    plot(kp_a,gmax,'DisplayName',CONAME_A{j_}); hold on;
   subplot(122)
    % plot(kp_a,wmax,'DisplayName',CONAME_A{j_}); hold on;
    plot(kp_a,c44,'DisplayName',CONAME_A{j_}); hold on;
end
   subplot(121)
legend('show'); grid on;
ylim([-10,10]);
xlabel('$k_\perp$'); ylabel('$\gamma_{max}$ from Eig(iCa)')
   subplot(122)
legend('show'); grid on;
xlabel('$k_\perp$'); ylabel('$\omega_{max}$ from Eig(iCa)')
ylim([-10,10]);

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
        '/home/ahoffman/gyacomo/iCa/gk_landau_P12_J6_dk_5e-2_km_1.0_NFLR_16.h5';...
        '/home/ahoffman/gyacomo/iCa/gk_landauii_P16_J9_dk_5e-2_km_2.0_NFLR_8.h5';...
        % '/home/ahoffman/gyacomo/iCa/gk_sugama_tau1e-3_P11_J7_dk_5e-2_km_5.0_NFLR_12.h5';...
        % '/home/ahoffman/gyacomo/iCa/gk_sugama_tau1e-3_P4_J2_dk_5e-2_km_5.0_NFLR_5.h5';...
        % '/home/ahoffman/gyacomo/iCa/gk_sugama_tau1e-2_P4_J2_dk_5e-2_km_5.0_NFLR_5.h5';...
        };
CONAME_A = {'SG 20 10';...
    'PA 20 10';...
    'FC 4 2 NFLR 6';...
    'FC 4 2 NFLR 12';...
    'LD 10 5 NFLR 12';...
    'LD 12 6 NFLR 12';...
    'LD 16 9 NFLR 16';...
    % 'SG 11 7 NFLR 12, tau 1e-3'; ...
    % 'SG 4 2 NFLR 5, tau 1e-3'; ...
    % 'SG 4 2 NFLR 5, tau 1e-2'; ...
    };
TAU_A = [1;...
    1;...
    1;...
    1;...
    1;...
    1;...
    % 1e-3;...
    % 1e-3;...
    % 1e-2;...
    ];
grow = {};
puls = {};
for j_ = 1:numel(mfns)
    mat_file_name = mfns{j_}; disp(mat_file_name);
    kp_a          =  h5read(mat_file_name,'/coordkperp');
    [~,idx_]       = min(abs(kp_a-kperp));
    matidx        = sprintf('%5.5i',idx_);
    FC_MAT           = h5read(mat_file_name,['/',matidx,'/Caapj/Ciipj']);
    grow{j_}  = real(eig(FC_MAT))*TAU_A(j_); 
    puls{j_}  = imag(eig(FC_MAT))*TAU_A(j_); 
end

figure
for j_ = 1:numel(mfns)
%    plot(puls{j_}, grow{j_},'o','DisplayName',CONAME_A{j_}); hold on;
   plot(grow{j_},'o','DisplayName',CONAME_A{j_}); hold on;
end
legend('show'); grid on; title(['$k_\perp=$',num2str(kperp)]);
xlabel('$\omega$ from Eig(iCa)'); ylabel('$\gamma$ from Eig(iCa)')
end