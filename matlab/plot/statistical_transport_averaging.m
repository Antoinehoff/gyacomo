function [ fig ] = statistical_transport_averaging( data, options )
scale = (1/data.Nx/data.Ny)^2;
Trange  = options.T;
[~,it0] = min(abs(Trange(1)-data.Ts0D)); 
[~,it1] = min(abs(Trange(end)-data.Ts0D)); 
gamma   = data.PGAMMA_RI(it0:it1)*scale;
dt_const = numel(unique(round(diff(data.Ts0D(it0:it1))*100)))==1;
if ~dt_const
    disp('DT not const on given interval');
else
    
    Ntot = (it1-it0)+1;

    transp_seg_avg = 1:Ntot;
    transp_seg_std = 1:Ntot;

    for Np = 1:Ntot % Loop on the number of segments
        transp_seg_avg(Np) = mean(gamma(1:Np));
        transp_seg_std(Np) = std(gamma(1:Np));
    end

    time_seg = (data.Ts0D(it0:it1)-data.Ts0D(it0)); 

    fig = figure;
%     subplot(211)
    plot(time_seg,transp_seg_avg,'-'); hold on;
    xlabel('Averaging time'); ylabel('transport value');
    
%     subplot(212)
%     errorbar(N_seg,transp_seg_avg,transp_seg_std);
%     xlabel('Averaging #points'); ylabel('transport value');
end
end

