DIR = '/misc/gyacomo23_outputs/scaling/strong/17x9x256x256x32/'; scaling='strong';
% DIR = '/misc/gyacomo23_outputs/scaling/strong/9x5x768x192x32/'; scaling='strong';
% DIR = '/misc/gyacomo23_outputs/scaling/strong/3x2x768x192x32/'; scaling='strong';
% DIR = '/misc/gyacomo23_outputs/scaling/weak/Np_5x2x128x64x32/'; scaling='weak';
% DIR = '/misc/gyacomo23_outputs/scaling/weak/Ny_3x2x768x8x32/'; scaling='weak';
% DIR = '/misc/gyacomo23_outputs/scaling/weak/Nz_3x2x768x32x8/'; scaling='weak';
% DIR = '/misc/gyacomo23_outputs/scaling/weak/Nz_5x2x128x64x8/'; scaling='weak';
% DIR = '/misc/gyacomo23_outputs/scaling/weak/Npyz_4x2x32x16x16/'; scaling='weak';

% Get a list of all items in the current directory
contents = dir(DIR);
Ncont = length(contents);

% Iterate through the contents to take valid dirs
dirs_ = {};
outp_ = {}; Noutp = 0;
for i = 1:Ncont
    subdir = [DIR,contents(i).name];
    filename = [subdir,'/outputs_00.h5'];
    dirs_{i} = subdir;
    try
        Np = h5readatt(filename,'/data/input/parallel','num_procs_p');
        Noutp  = Noutp + 1;
        outp_{Noutp} = filename;
    catch
    end
end

Np     = 0*(1:Noutp);
Ny     = 0*(1:Noutp);
Nz     = 0*(1:Noutp);
Nvar   = 0*(1:Noutp);
Rt_avg = 0*(1:Noutp);
Rt_std = 0*(1:Noutp);

% Iterate through the contents
for i = 1:Noutp
    % Get and display the name of the subdirectory
    filename = outp_{i};
    Np(i)  = h5readatt(filename,'/data/input/parallel','num_procs_p');
    Ny(i)  = h5readatt(filename,'/data/input/parallel','num_procs_ky');
    Nz(i)  = h5readatt(filename,'/data/input/parallel','num_procs_z');

    inp    = read_namelist([dirs_{i},'/fort_00.90']);
    Nvar(i)= (inp.GRID.pmax+1)*(inp.GRID.jmax+1)*inp.GRID.Nx*inp.GRID.Ny*inp.GRID.Nz;

    CPUTI = double(h5readatt(filename,'/data/input','cpu_time'));
    DTSIM = h5readatt(filename,'/data/input/basic','dt');

    RHSTC = h5read(filename,'/profiler/Tc_rhs');
    POITC = h5read(filename,'/profiler/Tc_poisson');
    SAPTC = h5read(filename,'/profiler/Tc_Sapj');
    EXBTC = h5read(filename,'/profiler/Tc_ExBshear');
    COLTC = h5read(filename,'/profiler/Tc_coll');
    GRATC = h5read(filename,'/profiler/Tc_grad');
    NADTC = h5read(filename,'/profiler/Tc_nadiab');
    ADVTC = h5read(filename,'/profiler/Tc_adv_field');
    GHOTC = h5read(filename,'/profiler/Tc_ghost');
    CLOTC = h5read(filename,'/profiler/Tc_clos');
    CHKTC = h5read(filename,'/profiler/Tc_checkfield');
    STETC = h5read(filename,'/profiler/Tc_step');
    TS0TC = h5read(filename,'/profiler/time');
    DIATC = h5read(filename,'/profiler/Tc_diag');
    CPUTI = CPUTI(end);

    RT_   = RHSTC + POITC + SAPTC + COLTC + GRATC + NADTC + ADVTC + ...
            GHOTC + CLOTC + CHKTC;
    if(numel(RT_)>5)
        Rt_avg(i) = mean(diff(RT_));
        Rt_std(i) = std(diff(RT_));
        % tmp_ = diff(RT_); Rt_avg(i) = tmp_(end);
    end
end
Np_tot = Np.*Ny.*Nz;
Rt_ref = Rt_avg(Np_tot==1);
Np_ref = 1;
MaxN  = max(Np_tot);
Np1Nz1 = logical((Np==1).*(Nz==1));
Np2Nz1 = logical((Np==2).*(Nz==1));
Np1Nz2 = logical((Np==1).*(Nz==2));
Np1Nz4 = logical((Np==1).*(Nz==4));
%%
figure;  hold on;
xlabel('Number of cores');
switch scaling
case 'strong'
    plot(Np_tot(Np1Nz1),Rt_ref./Rt_avg(Np1Nz1),...
        'o','DisplayName','$N_{pp}=1, N_{pz}=1$');
    plot(Np_tot(Np2Nz1),Rt_ref./Rt_avg(Np2Nz1), ...
        'o','DisplayName','$N_{pp}=2, N_{pz}=1$');
    plot(Np_tot(Np1Nz2),Rt_ref./Rt_avg(Np1Nz2), ...
        'o','DisplayName','$N_{pp}=1, N_{pz}=2$');
    plot(Np_tot(Np1Nz4),Rt_ref./Rt_avg(Np1Nz4), ...
        'o','DisplayName','$N_{pp}=1, N_{pz}=4$');
    plot(1:MaxN,(1:MaxN)/Np_ref,'--k','DisplayName','Ideal');
    set(gca,'XScale','log'); set(gca,'YScale','log');
    ylabel('Speedup');
    % title('Strong scaling speedup')
case 'weak'
    hold on;
    filt = logical((Np==2).*(Ny>=1).*(Nz>=1));
    plot(Np_tot(filt),Rt_ref./Rt_avg(filt),...
        'o','DisplayName','Effective');
    plot(1:MaxN,ones(1,MaxN),'--k','DisplayName','Ideal');
    ylim([0 1]);
    set(gca,'XScale','log')
     ylabel('Efficiency');
    % title('Weak scaling efficiency')
end