%% HeLaZ data
filename = 'results_00.h5';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data
moment = 'Ni00';

kr       = h5read(filename,['/data/var2d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var2d/' moment '/coordkz']);
timeNi   = h5read(filename,'/data/var2d/time');
Nipj     = zeros(numel(timeNi),numel(kr),numel(kz));
for it = 1:numel(timeNi)
    tmp          = h5read(filename,['/data/var2d/', moment,'/', num2str(it,'%06d')]);
    Nipj(it,:,:) = tmp.real + 1i * tmp.imaginary; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot growth rate vs kz
K_RICCI = 1; %% add a sqrt(1+tau) to the kperps
ikr     = 1; %% Fix the kr value
plot_gamma_vs_k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Ni00 evolution
ikr     = 1; %% Fix the kr value
plot_Ni00_t_evolution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%