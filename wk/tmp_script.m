% CO_A = {'DG','SG'};
% % CO_A = {'LD'};
% P_A  = [4 8 16 32];
% 
% for ic = 1:numel(CO_A)
%     for ip = 1:numel(P_A)
%         CO = CO_A{ic};
%         P  = P_A(ip);
%         CBC_kT_nu_scan
%     end
% end

%% Load and analyse
P_A  = [4 8 16 32];
% CO = 'DG';
CO = 'SG';
for i = 1:numel(P_A)
datafname = ['p2_linear_new/8x24_ky_0.3_P_',num2str(P_A(i)),...
    '_J_',num2str(P_A(i)/2),'_kT_2.5_5_nu_0.001_1_',CO,'GK.mat'];
load_metadata_scan
end