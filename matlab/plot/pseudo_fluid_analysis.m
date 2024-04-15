%% Pseudo fluid analysis
[Napjz, T_] = compile_results_3Da(DATADIR,J0,J1,'Napjz');
[~,it1] = min(abs(T_-Time_window(1)));
[~,it2] = min(abs(T_-Time_window(2)));

Napjz  = squeeze(Napjz(:,:,:,it1:it2));
E_dens = sum(abs(Napjz(1,1,1,:)));
E_Upar = sum(abs(Napjz(1,2,1,:)));
% E_Tpar = sum(abs(Napjz(1,3,1,:)));
E_Tpar = sum(abs(sqrt(2)*Napjz(1,3,1,:) + Napjz(1,3,1,:)));
E_Tper = sum(abs(Napjz(1,1,2,:)));
E_qpar = sum(abs(Napjz(1,4,1,:)));
E_qper = sum(abs(Napjz(1,2,2,:)));

E_i = [E_dens; E_Upar; E_Tpar; E_Tper; E_qpar; E_qper];

E_dens = sum(abs(Napjz(2,1,1,:)));
E_Upar = sum(abs(Napjz(2,2,1,:)));
E_Tpar = sum(abs(Napjz(2,3,1,:)));
E_Tper = sum(abs(Napjz(2,1,2,:)));
E_qpar = sum(abs(Napjz(2,4,1,:)));
E_qper = sum(abs(Napjz(2,2,2,:)));

E_e = [E_dens; E_Upar; E_Tpar; E_Tper; E_qpar; E_qper];

% E_tot = sum(E_e)+sum(E_i);
E_tot = sum(abs(Napjz(:)));
E_i = E_i/E_tot*100;
E_e = E_e/E_tot*100;

nms = {'$N_a^{00}$'; '$N_a^{10}$'; '$N_a^{20}$'; '$N_a^{01}$';...
    '$N_a^{30}$'; '$N_a^{11}$'};

figure;
bar([E_e E_i],'stacked'); hold on
% pie(E_i,nms); hold on
xticklabels(nms);
ylabel('$E [\%]$')
title(['$(',num2str(data.inputs.PMAX),',',num2str(data.inputs.JMAX),')$',...
    ', $k_y=$',num2str(data.grids.ky(2)),...
    ', $s=',num2str(sum(E_i+E_e)),'\%$']);
