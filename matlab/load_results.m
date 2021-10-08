%% load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Loading ',filename])
% Loading from output file
CPUTIME   = h5readatt(filename,'/data/input','cpu_time');
DT_SIM    = h5readatt(filename,'/data/input','dt');

[Pe, Je, Pi, Ji, kx, ky, z] = load_grid_data(filename);

W_GAMMA   = strcmp(h5readatt(filename,'/data/input','write_gamma'),'y');
W_HF      = strcmp(h5readatt(filename,'/data/input','write_hf'   ),'y');
W_PHI     = strcmp(h5readatt(filename,'/data/input','write_phi'  ),'y');
W_NA00    = strcmp(h5readatt(filename,'/data/input','write_Na00' ),'y');
W_NAPJ    = strcmp(h5readatt(filename,'/data/input','write_Napj' ),'y');
W_SAPJ    = strcmp(h5readatt(filename,'/data/input','write_Sapj' ),'y');
W_DENS    = strcmp(h5readatt(filename,'/data/input','write_dens' ),'y');
W_TEMP    = strcmp(h5readatt(filename,'/data/input','write_temp' ),'y');


if W_GAMMA
    [ GGAMMA_RI, Ts0D, dt0D] = load_0D_data(filename, 'gflux_ri');
      PGAMMA_RI              = load_0D_data(filename, 'pflux_ri');
end

if W_HF
    [ HFLUX_X, Ts0D, dt0D] = load_0D_data(filename, 'hflux_x');
end

if W_PHI
    [ PHI, Ts3D, dt3D] = load_3D_data(filename, 'phi');
end

if W_NA00
    [Ni00, Ts3D, dt3D] = load_3D_data(filename, 'Ni00');
     Ne00              = load_3D_data(filename, 'Ne00');
end


if W_NAPJ
    [Nipj, Ts5D, dt5D] = load_5D_data(filename, 'moments_i');
    [Nepj            ] = load_5D_data(filename, 'moments_e');
end

if W_SAPJ
    [Sipj, Ts5D, dt5D] = load_5D_data(filename, 'Sipj');
     Sepj              = load_5D_data(filename, 'Sepj');
end

if W_DENS
    [DENS_E, Ts3D, dt3D] = load_3D_data(filename, 'dens_e');
    [DENS_I, Ts3D, dt3D] = load_3D_data(filename, 'dens_i');
end

if W_TEMP
    [TEMP_E, Ts3D, dt3D] = load_3D_data(filename, 'temp_e');
    [TEMP_I, Ts3D, dt3D] = load_3D_data(filename, 'temp_i');
end