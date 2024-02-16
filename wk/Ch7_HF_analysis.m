PARTITION = '/misc/gyacomo23_outputs/paper_3/';
switch 3
case 1
    SIM_SET_NAME = 'Multi-scale';
    E_FLUX       = 1;
    FOLDER       = 'DIIID_fullphys_rho95_geom_scan/multi_scale/multi_scale_3x2x768x192x24';
    % FOLDER       = 'DIIID_fullphys_rho95_geom_scan/multi_scale/multi_scale_5x2x768x192x24';
case 2
    SIM_SET_NAME = 'Ion-scale';
    E_FLUX       = 1;
    FOLDER       = 'DIIID_fullphys_rho95_geom_scan/ion_scale/ion_scale_5x2x256x64x32_tau_1_RN_0';
case 3
    SIM_SET_NAME = 'Adiab. e.';
    E_FLUX       = 0;
    FOLDER       = 'DIIID_adiab_e_rho95_geom_scan/5x2x256x64x32_tau_1_RN_0/';
case 4
    SIM_SET_NAME = 'Adiab. e.';
    E_FLUX       = 0;
    FOLDER       = 'DIIID_adiab_e_rho95_geom_scan/5x2x256x64x32_tau_1_RN_0/';    
end

GEOM = {'NT','0T','PT'};
COLO = {'b','k','r'};
figure
for i = 1:3
    OPTION = GEOM{i};
    datadir = [PARTITION,FOLDER,'/',OPTION];
    out = read_flux_out_XX(datadir,0);
    
    [ts, is] = sort(out.t);
    Pxis = out.Pxi(is);
    Qxis = out.Qxi(is);

    if E_FLUX
        Pxes = out.Pxe(is);
        Qxes = out.Qxe(is);
        subplot(211)
        plot(ts,Qxes,COLO{i}); hold on;
        subplot(212)
        title(SIM_SET_NAME)
    end
    plot(ts,Qxis,COLO{i}); hold on;
end
title(SIM_SET_NAME)

