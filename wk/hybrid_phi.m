%% Load HeLaZ results
% Nepj

moment = 'moments_e';

kr       = h5read(filename,['/data/var5d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var5d/' moment '/coordkz']);
time     = h5read(filename,'/data/var5d/time');
Nepj     = zeros(GRID.pmaxe+1, GRID.jmaxe+1,numel(kr),numel(kz),numel(time));
for it = 1:numel(timeNe)
    tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
    Nepj(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary; 
end

% Nipj
moment = 'moments_i';

timeNi   = h5read(filename,'/data/var5d/time');
Nipj     = zeros(GRID.pmaxi+1, GRID.jmaxi+1,numel(kr),numel(kz),numel(time));
for it = 1:numel(timeNe)
    tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
    Nipj(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary; 
end

NapjHeLaZ = zeros(size(results.Napj));
NapjHeLaZ(:, 1:bare(GRID.pmaxe, GRID.jmaxe))     = transpose(reshape(Nepj,[numel(time),4]));
NapjHeLaZ(:, bare(GRID.pmaxe, GRID.jmaxe)+1:end) = transpose(reshape(Nipj,[numel(time),4]));
%% MOLI results
NepjMOLI = results.Napj(:, 1:bare(GRID.pmaxe, GRID.jmaxe)    );
NipjMOLI = results.Napj(:, bare(GRID.pmaxe, GRID.jmaxe)+1:end);

%% phi
time  = h5read(filename,'/data/var2d/time');
kr       = h5read(filename,'/data/var2d/phi/coordkr');
kz       = h5read(filename,'/data/var2d/phi/coordkz');
phiHeLaZ      = zeros(numel(time),numel(kr),numel(kz));
for it = 1:numel(time)
    tmp         = h5read(filename,['/data/var2d/phi/' num2str(it,'%06d')]);
    phiHeLaZ(it,:,:) = tmp.real + 1i * tmp.imaginary;
end

% Phi from MOLI routine with HeLaZ results
phihybrid = zeros(size(time));
for it = 1:numel(time)
    phihybrid(it) = get_phi(NapjHeLaZ(it,:),params,options);
end

% Phi MOLI
phiMOLI  = zeros(size(time));
for it = 1:numel(time)
    phiMOLI(it) = get_phi(results.Napj(it,:),params,options);
end

fig = figure;
%HeLaZ results
semilogy(time,abs(phiHeLaZ),'-','DisplayName','HeLaZ RK4')
hold on
title(TITLE);
%hybrid results
semilogy(time,abs(phihybrid),'-.','DisplayName','hybrid ')
%MOLI results
semilogy(time,abs(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
grid on
xlabel('$t$')
ylabel('$|\phi|$')
legend('show')
if SAVEFIG; FIGNAME = ['phi_kz_',num2str(GRID.kzmin)]; save_figure; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%