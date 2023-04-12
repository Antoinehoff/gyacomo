clrs_ = lines(10);
kN=2.22;
figure
for i = 1:3
if i==1
prefix = '/home/ahoffman/gyacomo/results/paper_2_GYAC23/3x2x64x48x16/';
resdirs = {...
    'CBC';...
    'kT_5.3';...
    'kT_4.5';...
    'kT_4.0';...
    'kT_3.5';...
    'kT_3.0';...
    };
pname = '3x2x64x48x16';  clr_ = clrs_(1,:);
ITG_threshold = 3.2;
end

if i==2
prefix = '/home/ahoffman/gyacomo/results/paper_2_GYAC23/5x2x64x48x16/';
resdirs = {...
    'CBC';...
    'kT_5.3';...
    'kT_4.5';...
    'kT_4.0';...
    'kT_3.5';...
    'kT_3.0';...
    };
pname = '5x2x64x48x16'; clr_ = clrs_(2,:);
ITG_threshold = 3.7;
end

if i==3
prefix = '/home/ahoffman/gyacomo/results/paper_2_GYAC23/9x2x64x48x16/';
resdirs = {...
    'CBC';...
    'kT_5.3';...
    'kT_4.5';...
    'kT_4.0';...
    'kT_3.5';...
    'kT_3.0';...
    };
pname = '9x2x64x48x16'; clr_ = clrs_(3,:);
ITG_threshold = 4.2;
end

if i==4
    
ITG_threshold > 3.6; 
end

x = [...
    6.96,...
    5.3,...
    4.5,...
    4.0,...
    3.5,...
    3.0...
    ];

J0 = 00; J1 = 10;

Nseg = 5;

Qx_avg = 0*(1:numel(resdirs));
Qx_std = 0*(1:numel(resdirs));

for i = 1:numel(resdirs)
    data    = compile_results_low_mem(data,[prefix,resdirs{i},'/'],J0,J1);
    Trange  = data.Ts0D(end)*[0.3 1.0];
    %
    [~,it0] = min(abs(Trange(1)  -data.Ts0D)); 
    [~,it1] = min(abs(Trange(end)-data.Ts0D)); 
    %
    if 0
        Qx      = data.HFLUX_X(it0:it1);
        Qxa_    = 0*(1:Nseg);
        for n = 1:Nseg
           ntseg = floor((it1-it0)/n);
           for m = 1:n 
            Qxa_(n) = Qxa_(n) + mean(Qx((1:ntseg)+(m-1)*ntseg));
           end
           Qxa_(n) = Qxa_(n)/n;
        end
        Qx_avg(i) = mean(Qxa_);
        Qx_std(i) = std(Qxa_);
    else
        Qx_avg(i) = mean(data.HFLUX_X(it0:it1));
        Qx_std(i) =  std(data.HFLUX_X(it0:it1));
    end
end
Chi_avg = Qx_avg./x/kN;
Chi_std = Qx_std./x/kN;
% plot;
errorbar(x,Chi_avg,Chi_std,'DisplayName',pname,'color',clr_); hold on;
plot(ITG_threshold*[1 1],[0 20],'-.','DisplayName','$\kappa_T^{crit}$',...
    'color',clr_);
end
ylabel('$\chi$');
xlabel('$\kappa_T (\kappa_N=2.22)$');
ylim([0,10]);
legend('show');
