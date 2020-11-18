%% load profiling
filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],00);

CPUTIME   = h5readatt(filename,'/data/input','cpu_time');
DT_SIM    = h5readatt(filename,'/data/input','dt');


rhs_Tc       = h5read(filename,'/profiler/Tc_rhs');
adv_field_Tc = h5read(filename,'/profiler/Tc_adv_field');
poisson_Tc   = h5read(filename,'/profiler/Tc_poisson');
Sapj_Tc      = h5read(filename,'/profiler/Tc_Sapj');
diag_Tc      = h5read(filename,'/profiler/Tc_diag');
checkfield_Tc= h5read(filename,'/profiler/Tc_checkfield');
step_Tc      = h5read(filename,'/profiler/Tc_step');
Ts0D         = h5read(filename,'/profiler/time');

missing_Tc   = step_Tc - rhs_Tc - adv_field_Tc -...
               poisson_Tc - Sapj_Tc -diag_Tc -checkfield_Tc;
total_Tc     = step_Tc;
%% Plots
Y = [diff(rhs_Tc); diff(adv_field_Tc); diff(poisson_Tc);...
    diff(Sapj_Tc); diff(checkfield_Tc); diff(diag_Tc); diff(missing_Tc)];
Y = reshape(Y,[numel(Y)/7,7]);
fig = figure;

p1 = area(Ts0D(2:end),100*Y./diff(total_Tc),'LineStyle','none');
legend('Compute RHS','Adv. fields','Poisson','Sapj','Check+Sym','Diag','Missing')
xlabel('Sim. Time'); ylabel('Step Comp. Time [\%]')
ylim([0,100]); xlim([Ts0D(2),Ts0D(end)]);
hold on
yyaxis right
p2 = plot(Ts0D(2:end),diff(total_Tc),'--r','LineWidth',1.0);
ylabel('Step Comp. Time [s]')
ylim([0,1.1*max(diff(total_Tc))])
set(gca,'ycolor','r') 
FIGNAME = 'profiler';
save_figure