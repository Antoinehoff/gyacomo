function [] = profiler(data)
%% load profiling
% filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],00);
% outfilename = ['/misc/HeLaZ_outputs',filename(3:end)];
CPUTI=[]; DTSIM=[]; RHSTC=[]; POITC=[]; SAPTC=[]; COLTC=[];
GRATC=[]; NADTC=[];  ADVTC=[]; GHOTC=[]; CLOTC=[]; CHKTC=[];
DIATC=[]; STETC=[]; TS0TC=[];
for i = 1:numel(data.outfilenames)
    outfilename = data.outfilenames{i};
    CPUTI = [ CPUTI; double(h5readatt(outfilename,'/data/input','cpu_time'))];
    DTSIM = [ DTSIM; h5readatt(outfilename,'/data/input/basic','dt')];

    RHSTC = [ RHSTC; h5read(outfilename,'/profiler/Tc_rhs')];
    POITC = [ POITC; h5read(outfilename,'/profiler/Tc_poisson')];
    SAPTC = [ SAPTC; h5read(outfilename,'/profiler/Tc_Sapj')];
    COLTC = [ COLTC; h5read(outfilename,'/profiler/Tc_coll')];
    GRATC = [ GRATC; h5read(outfilename,'/profiler/Tc_grad')];
    NADTC = [ NADTC; h5read(outfilename,'/profiler/Tc_nadiab')];
    ADVTC = [ ADVTC; h5read(outfilename,'/profiler/Tc_adv_field')];
    GHOTC = [ GHOTC; h5read(outfilename,'/profiler/Tc_ghost')];
    CLOTC = [ CLOTC; h5read(outfilename,'/profiler/Tc_clos')];
    CHKTC = [ CHKTC; h5read(outfilename,'/profiler/Tc_checkfield')];
    DIATC = [ DIATC; h5read(outfilename,'/profiler/Tc_diag')];
    STETC = [ STETC; h5read(outfilename,'/profiler/Tc_step')];
    TS0TC = [ TS0TC; h5read(outfilename,'/profiler/time')];
end

N_T          = 12;

missing_Tc   = STETC - RHSTC - ADVTC - GHOTC - CLOTC ...
              -COLTC - POITC - SAPTC - CHKTC - DIATC - GRATC - NADTC;
total_Tc     = STETC;

TIME_PER_FCT = [diff(RHSTC); diff(ADVTC); diff(GHOTC);...
    diff(CLOTC); diff(COLTC); diff(POITC); diff(SAPTC); ...
    diff(CHKTC); diff(DIATC); diff(GRATC); diff(NADTC); diff(missing_Tc)];
TIME_PER_FCT = reshape(TIME_PER_FCT,[numel(TIME_PER_FCT)/N_T,N_T]);

TIME_PER_STEP = sum(TIME_PER_FCT,2);
TIME_PER_CPU  = trapz(TS0TC(2:end),TIME_PER_STEP);

rhs_Ta        = mean(diff(RHSTC));
adv_field_Ta  = mean(diff(ADVTC));
ghost_Ta      = mean(diff(GHOTC));
clos_Ta       = mean(diff(CLOTC));
coll_Ta       = mean(diff(COLTC));
poisson_Ta    = mean(diff(POITC));
Sapj_Ta       = mean(diff(SAPTC));
checkfield_Ta = mean(diff(CHKTC));
grad_Ta       = mean(diff(GRATC));
nadiab_Ta     = mean(diff(NADTC));
diag_Ta       = mean(diff(DIATC));
miss_Ta       = mean(diff(missing_Tc));
total_Ta      = mean(diff(total_Tc));
names = {...
    'Mrhs';
    'Advf';
    'Ghst';
    'Clos';
    'Capj';
    'Pois';
    'Sapj';
    'Chck';
    'Diag';
    'Grad';
    'napj';
    'Miss';
};
Ts_A = [RHSTC(end) ADVTC(end) GHOTC(end) CLOTC(end) COLTC(end)...
    POITC(end) SAPTC(end) CHKTC(end) DIATC(end) GRATC(end) NADTC(end) missing_Tc(end)];
NSTEP_PER_SAMP= mean(diff(TS0TC))/DTSIM;

%% Plots
if 1
    %% Area plot
fig = figure;
% colors = rand(N_T,3);
% colors = lines(N_T);
colors = distinguishable_colors(N_T);
x_ = TS0TC(2:end);
y_ = TIME_PER_FCT;
xx_= zeros(2*numel(x_),1);
yy_= zeros(2*numel(x_),numel(names));
dx = (x_(2)-x_(1));
xx_(1) = x_(1)-dx/2; xx_(2) = x_(1)+dx/2;
yy_(1,:) = y_(1,:)/dx;     
yy_(2,:) = y_(2,:)/dx;
for i = 2:numel(x_)
    dx = x_(i) - x_(i-1);
    xx_(2*i-1) = x_(i)-dx/2;
    xx_(2*i  ) = x_(i)+dx/2;
    yy_(2*i-1,:) = y_(i,:)/dx;
    yy_(2*i  ,:) = y_(i,:)/dx;
end
p1 = area(xx_,yy_,'LineStyle','none');
for i = 1:N_T; p1(i).FaceColor = colors(i,:);
    LEGEND{i} = sprintf('%s t=%1.1e[s] (%0.1f %s)',names{i},Ts_A(i),Ts_A(i)/total_Tc(end)*100,'\%');
end;
legend(LEGEND,'Location','bestoutside');
% legend('Compute RHS','Adv. fields','ghosts comm', 'closure', 'collision','Poisson','Nonlin','Check+sym', 'Diagnos.', 'Process', 'Missing')
xlabel('Sim. Time [$\rho_s/c_s$]'); ylabel('N Steps Comp. Time [s]')
xlim([TS0TC(2),TS0TC(end)]);
title(sprintf('Gyacomo 2, CPU time:  ~%.0f [h] ~%.0f [min] ~%.0f [s]',...
    CPUTI/3600,CPUTI/60-60*floor(CPUTI/3600),CPUTI-60*floor(CPUTI/60-60*floor(CPUTI/3600))))
hold on
FIGNAME = 'profiler';
% save_figure
% Re-order Legend
lbl = fig.Children(1).String;                         % Retrieve legend labels
numlbl = length(lbl);                               % Determine number of lables
order = sort(1:1:numlbl,'descend');                 % Create array of label numbers in descending order
newlbl = lbl(order);                                % Create new labels in descending order
legend(findobj(fig.Children(2),'Type','area'),newlbl) % Set the legend to follow the new labels
else
    %% Normalized Area plot
fig = figure;
colors = colorcube(N_T);
p1 = area(Ts0D(2:end),100*TIME_PER_FCT./diff(total_Tc),'LineStyle','none', 'FaceColor','flat');
% for i = 1:N_T; p1(i).FaceColor = rand(1,3); end;
for i = 1:N_T; p1(i).FaceColor = colors(i,:); end;
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
% save_figure
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
histogram(diff(process_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(checkfield_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(diag_Tc)/NSTEP_PER_SAMP,'Normalization','probability');hold on
grid on;
legend('Compute RHS','Adv. fields','Ghosts comm', 'closure', 'collision','Poisson','Nonlin','Process','Check+sym', 'Diagnos.', 'Missing')
xlabel('Step Comp. Time [s]'); ylabel('')
set(gca,'Xscale','log')
FIGNAME = 'profiler';
% save_figure
end
end