
% Nominal parameters
PARTITION = '/misc/gyacomo23_outputs/triangularity_paper/';
% resdir = 'ion_scale/3x2x256x64x32/0T';
% resdir = 'ion_scale/5x2x256x64x32/NT';
% resdir = 'ion_scale/5x3x256x64x32/NT';
% resdir = 'ion_scale/5x3x192x48x24/NT';
resdir = 'ion_scale/9x5x256x64x32/0T';
% resdir = 'ion_scale/restart/5x3x256x64x32/0T';
% resdir = 'ion_scale/restart/9x5x192x48x24';
% resdir = 'adiabatic_electrons/5x2x256x64x32/0T';
% resdir = 'adiabatic_electrons/5x2x192x48x24/0T';
% resdir = 'hot_electrons/256x64x32/0T';
% resdir = 'hot_electrons/256x64x32/0T';
% resdir = 'hot_electrons/512x64x32/0T';
DATADIR = [PARTITION,resdir,'/'];
read_flux_out_XX(DATADIR,1,1);
%%
% No grad N
PARTITION = '/misc/gyacomo23_outputs/triangularity_paper/no_gradN/';
% resdir = '/ion_scale/3x2x256x64x24/0T';
% resdir = '/ion_scale/5x2x256x64x24/PT';
% resdir = '/ion_scale/5x3x256x64x24/PT';
% resdir = '/ion_scale/5x3x192x48x24/NT';
% resdir = '/adiabatic_electrons/3x2x256x64x24/0T';
% resdir = '/adiabatic_electrons/5x3x192x48x24/0T';
% resdir = '/adiabatic_electrons/5x2x256x64x24/NT';
resdir = '/hot_electrons/256x64x32/0T/noise_init';
% resdir = '/hot_electrons/L_300/256x64x32/NT';
DATADIR = [PARTITION,resdir,'/'];
read_flux_out_XX(DATADIR,1,10);
%%
PARTITION = '/misc/gyacomo23_outputs/triangularity_paper/no_gradN/';
Models = {'ion_scale/5x3x192x48x24',...
    'adiabatic_electrons/5x3x192x48x24',...
    'hot_electrons/256x64x32'};
% Models = {'ion_scale/5x2x256x64x24',...
    % 'adiabatic_electrons/5x2x256x64x24',...
    % 'hot_electrons/256x64x32'};
clrs   = lines(3);
triangularities = {'NT','0T','PT'};

figure
t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:3
    nexttile
    for j = 1:3
        DATADIR = [PARTITION,Models{j},'/',triangularities{i},'/'];
        out = read_flux_out_XX(DATADIR,0,1);
        plot(out.t, out.Qxi,'color',clrs(j,:)); hold on;
    end
    ylim([0 200]); ylabel('$Q_{xi}$');
    xlim([0 300]); xlabel('$tc_s/R$');
end

