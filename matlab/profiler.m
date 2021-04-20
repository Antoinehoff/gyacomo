%% load profiling
% filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],00);

CPUTIME   = double(h5readatt(filename,'/data/input','cpu_time'));
DT_SIM    = h5readatt(filename,'/data/input','dt');


rhs_Tc       = h5read(filename,'/profiler/Tc_rhs');
adv_field_Tc = h5read(filename,'/profiler/Tc_adv_field');
ghost_Tc      = h5read(filename,'/profiler/Tc_ghost');
clos_Tc      = h5read(filename,'/profiler/Tc_clos');
coll_Tc      = h5read(filename,'/profiler/Tc_coll');
poisson_Tc   = h5read(filename,'/profiler/Tc_poisson');
Sapj_Tc      = h5read(filename,'/profiler/Tc_Sapj');
checkfield_Tc= h5read(filename,'/profiler/Tc_checkfield');
diag_Tc      = h5read(filename,'/profiler/Tc_diag');
step_Tc      = h5read(filename,'/profiler/Tc_step');
Ts0D         = h5read(filename,'/profiler/time');

missing_Tc   = step_Tc - rhs_Tc - adv_field_Tc - ghost_Tc -clos_Tc ...
              -coll_Tc -poisson_Tc -Sapj_Tc -checkfield_Tc -diag_Tc;
total_Tc     = step_Tc;

TIME_PER_FCT = [diff(rhs_Tc); diff(adv_field_Tc); diff(ghost_Tc);...
    diff(clos_Tc); diff(coll_Tc); diff(poisson_Tc); diff(Sapj_Tc); ...
    diff(checkfield_Tc); diff(diag_Tc); diff(missing_Tc)];
TIME_PER_FCT = reshape(TIME_PER_FCT,[numel(TIME_PER_FCT)/10,10]);

TIME_PER_STEP = sum(TIME_PER_FCT,2);
TIME_PER_CPU  = trapz(Ts0D(2:end),TIME_PER_STEP);

rhs_Ta        = mean(diff(rhs_Tc));
adv_field_Ta  = mean(diff(adv_field_Tc));
ghost_Ta      = mean(diff(ghost_Tc));
clos_Ta       = mean(diff(clos_Tc));
coll_Ta       = mean(diff(coll_Tc));
poisson_Ta    = mean(diff(poisson_Tc));
Sapj_Ta       = mean(diff(Sapj_Tc));
checkfield_Ta = mean(diff(checkfield_Tc));
diag_Ta       = mean(diff(diag_Tc));

NSTEP_PER_SAMP= mean(diff(Ts0D))/DT_SIM;

%% Plots
if 1
    %% Area plot
fig = figure;

p1 = area(Ts0D(2:end),TIME_PER_FCT,'LineStyle','none');
for i = 1:10; p1(i).FaceColor = rand(1,3); end;
legend('Compute RHS','Adv. fields','ghosts comm', 'closure', 'collision','Poisson','Nonlin','Check+sym', 'Diagnos.', 'Missing')
xlabel('Sim. Time [$\rho_s/c_s$]'); ylabel('Step Comp. Time [s]')
xlim([Ts0D(2),Ts0D(end)]);
title(sprintf('Proc. 1, total sim. time  ~%.0f [h]',CPUTIME/3600))
hold on
FIGNAME = 'profiler';
save_figure

else
    %% Normalized Area plot
fig = figure;

p1 = area(Ts0D(2:end),100*TIME_PER_FCT./diff(total_Tc),'LineStyle','none', 'FaceColor','flat');
for i = 1:10; p1(i).FaceColor = rand(1,3); end;
legend('Compute RHS','Adv. fields','ghosts comm', 'closure', 'collision','Poisson','Nonlin','Check+sym', 'Diagnos.', 'Missing')
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
end

if 0
    %% Histograms
fig = figure;
histogram(diff(rhs_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(adv_field_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(ghost_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(clos_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(coll_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(poisson_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(Sapj_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(checkfield_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(diag_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
grid on;
legend('Compute RHS','Adv. fields','Ghosts comm', 'closure', 'collision','Poisson','Nonlin','Check+sym', 'Diagnos.', 'Missing')
xlabel('Step Comp. Time [s]'); ylabel('')
FIGNAME = 'profiler';
save_figure
end