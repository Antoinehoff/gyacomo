figure

PJ  = '5x3';
% PJ  = '9x5';
% PJ  = '11x6';
% PJ = '17x9';

% RES = '128x64x24';
RES = '192x96x24';

NU  = '0.005';

PJxRES = [PJ,'x',RES];
for CO = {'DGGK','SGGK','LDGK'}
folder = ['/misc/gyacomo23_outputs/paper_2_GYAC23/collision_study/nu',CO{1},'_scan_kT_5.3/',PJxRES,'/nu_',NU];
[t, ~, Qxi] = read_flux_out_XX(folder);
plot(t,Qxi,'.','DisplayName',CO{1}); hold on

end
legend('show');
title([PJ,'x',RES,' nu=',NU]);
ylim([0 25])