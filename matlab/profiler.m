function [] = profiler(DATADIR,JID)
%% load profiling
% filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],00);
% outfilename = ['/misc/HeLaZ_outputs',filename(3:end)];
CPUTI=[]; DTSIM=[]; RHSTC=[]; POITC=[]; SAPTC=[]; COLTC=[];
GRATC=[]; NADTC=[]; ADVTC=[]; GHOTC=[]; CLOTC=[]; CHKTC=[];
DIATC=[]; EXBTC=[]; STETC=[]; TS0TC=[];
% for i = 1:numel(data.outfilenames)
    outfilename = [DATADIR,'/',sprintf('outputs_%02d.h5',JID)];
    CPUTI = [ CPUTI; double(h5readatt(outfilename,'/data/input','cpu_time'))];
    DTSIM = [ DTSIM; h5readatt(outfilename,'/data/input/basic','dt')];

    RHSTC = [ RHSTC; h5read(outfilename,'/profiler/Tc_rhs')];
    POITC = [ POITC; h5read(outfilename,'/profiler/Tc_poisson')];
    SAPTC = [ SAPTC; h5read(outfilename,'/profiler/Tc_Sapj')];
    EXBTC = [ EXBTC; h5read(outfilename,'/profiler/Tc_ExBshear')];
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
% end
CPUTI = sum(CPUTI);
DTSIM = mean(DTSIM);
N_T          = 13;

MISTC   = STETC - RHSTC - ADVTC - GHOTC - CLOTC - EXBTC ...
              -COLTC - POITC - SAPTC - CHKTC - DIATC - GRATC - NADTC;
total_Tc     = STETC;

TIME_PER_FCT = [diff(RHSTC); diff(ADVTC); diff(GHOTC);...
    diff(CLOTC); diff(COLTC); diff(POITC); diff(SAPTC); diff(EXBTC); ...
    diff(CHKTC); diff(DIATC); diff(GRATC); diff(NADTC); diff(MISTC)];
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
ExB_Ta        = mean(diff(EXBTC));
checkfield_Ta = mean(diff(CHKTC));
grad_Ta       = mean(diff(GRATC));
nadiab_Ta     = mean(diff(NADTC));
diag_Ta       = mean(diff(DIATC));
miss_Ta       = mean(diff(MISTC));
total_Ta      = mean(diff(total_Tc));
names = {...
    'Mrhs';
    'Advf';
    'Ghst';
    'Clos';
    'Capj';
    'Pois';
    'Sapj';
    'ExBs';
    'Chck';
    'Diag';
    'Grad';
    'napj';
    'Miss';
};
Ts_A = [rhs_Ta adv_field_Ta ghost_Ta clos_Ta coll_Ta poisson_Ta...
     Sapj_Ta ExB_Ta checkfield_Ta diag_Ta grad_Ta  nadiab_Ta miss_Ta];
NSTEP_PER_SAMP= mean(diff(TS0TC))/DTSIM;

%% Plots
fig = figure;
if 1
    % subplot(121)
    %% Area plot
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
    yy_(2*i-1,:) = y_(i,:)/(dx/DTSIM);
    yy_(2*i  ,:) = y_(i,:)/(dx/DTSIM);
end
p1 = area(xx_/DTSIM,yy_,'LineStyle','none');
for i = 1:N_T; p1(i).FaceColor = colors(i,:);
%     LEGEND{i} = sprintf('%s t=%1.1e[s] (%0.1f %s)',names{i},Ts_A(i),Ts_A(i)/total_Tc(end)*100,'\%');
    LEGEND{i} = [names{i},' $\hat t=$',sprintf('%1.1e[s] (%0.1f %s)',Ts_A(i)/NSTEP_PER_SAMP,Ts_A(i)/total_Ta*100,'\%')];
end
legend(LEGEND,'Location','bestoutside');
% legend('Compute RHS','Adv. fields','ghosts comm', 'closure', 'collision','Poisson','Nonlin','Check+sym', 'Diagnos.', 'Process', 'Missing')
xlabel('Simulation step'); ylabel('Comp. time per step [s]')
xlim([TS0TC(2),TS0TC(end)]./DTSIM);
ylim([0, 1.5*total_Ta/(dx/DTSIM)])
h_ = floor(CPUTI/3600);
m_ = floor(floor(CPUTI/60)-60*h_);
s_ = CPUTI - 3600*h_ - 60*m_;
title([DATADIR,sprintf(' (%.0f [h] ~%.0f [min] ~%.0f [s])',...
    h_,m_,s_)])
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
    subplot(122)
    %% Histograms
% fig = figure;
histogram(diff(RHSTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(ADVTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(GHOTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(CLOTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(COLTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(POITC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(SAPTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(CHKTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(DIATC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
histogram(diff(MISTC)/NSTEP_PER_SAMP,'Normalization','probability');hold on
grid on;
legend('Compute RHS','Adv. fields','Ghosts comm', 'closure', 'collision','Poisson','Nonlin','Check+sym', 'Diagnos.', 'Missing')
xlabel('Step Comp. Time [s]'); ylabel('')
set(gca,'Xscale','log')
FIGNAME = 'profiler';
% save_figure
end
end