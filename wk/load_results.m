%% load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JOBNUM = 00;
FMT = '.fig';
filename = [BASIC.SIMID,'_','%.2d.h5'];
filename = sprintf(filename,JOBNUM); disp(['Loading ',filename])
% Loading from output file
% CPUTIME   = h5readatt(filename,'/data/input','cpu_time');
if strcmp(OUTPUTS.write_moments,'.true.') 
[Nipj, Pi, Ji, kr, kz, Ts5D, dt5D] = load_5D_data(filename, 'moments_i');
[Nepj, Pe, Je                    ] = load_5D_data(filename, 'moments_e');
end
[Ni00, kr, kz, Ts2D, dt2D] = load_2D_data(filename, 'Ni00');
 Ne00                      = load_2D_data(filename, 'Ne00');
%Pi = [0]; Ji = Pi; Pe = Pi; Je = Pi;
% Nipj = zeros(1,1,numel(kr),numel(kz),numel(Ts5D));
% Nepj = Nipj;
% Nipj(1,1,:,:,:) = Ni00; Nepj(1,1,:,:,:) = Ne00;
PHI                            = load_2D_data(filename, 'phi');

if strcmp(OUTPUTS.write_non_lin,'.true.') 
    Sipj    = load_5D_data(filename, 'Sipj');
    Sepj    = load_5D_data(filename, 'Sepj');
end