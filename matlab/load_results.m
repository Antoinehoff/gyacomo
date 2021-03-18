%% load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],JOBNUM);
disp(['Loading ',filename])
% Loading from output file
CPUTIME   = h5readatt(filename,'/data/input','cpu_time');
DT_SIM    = h5readatt(filename,'/data/input','dt');

[Pe, Je, Pi, Ji, kr, kz] = load_grid_data(filename);

if W_GAMMA
    [ GGAMMA_RI, Ts0D, dt0D] = load_0D_data(filename, 'gflux_ri');
      PGAMMA_RI              = load_0D_data(filename, 'pflux_ri');
end

if W_PHI
    [ PHI, Ts2D, dt2D] = load_2D_data(filename, 'phi');
end

if W_NA00
    [Ni00, Ts2D, dt2D] = load_2D_data(filename, 'Ni00');
     Ne00              = load_2D_data(filename, 'Ne00');
end


if W_NAPJ
    [Nipj, Ts5D, dt5D] = load_5D_data(filename, 'moments_i');
    [Nepj            ] = load_5D_data(filename, 'moments_e');
end

if W_SAPJ
    [Sipj, Ts5D, dt5D] = load_5D_data(filename, 'Sipj');
     Sepj              = load_5D_data(filename, 'Sepj');
end
