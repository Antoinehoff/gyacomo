gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options

PARTITION = '/Users/ahoffmann/gyacomo/results/triangularity_paper/no_gradN/';

models = {'ion_scale','adiabatic_electrons','hot_electrons'};
% resolu = {'5x3x192x48x24','5x3x192x48x24','256x64x32'};
% resolu = {'3x2x256x64x24','3x2x256x64x24','256x64x32'};
resolu = {'5x2x256x64x24','5x2x256x64x24','256x64x32'};
triang = {'NT','0T','PT'};
colors = lines(3);
figure
for i = 1:3
    for j = 1:3
        DATADIR = [PARTITION,'/',models{j},'/',resolu{j},'/',triang{i},'/'];
        out = read_flux_out_XX(DATADIR,0,10);
        subplot(1,3,i);
        plot(out.t,out.Qxi,'color',colors(j,:)); hold on;
        ylim([0 200]); xlim([0 300]);
    end
end
