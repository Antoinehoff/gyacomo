%% Default values for plots
set(0,'defaulttextInterpreter','latex') 
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)
set(0, 'DefaultLineLineWidth', 2.0);
line_colors = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];...
    [0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840]];
line_styles = {'-','-.','--',':'};

%Title management
if     MODEL.CO == -1; CONAME = 'FC';
elseif MODEL.CO == -2; CONAME = 'DC';
elseif MODEL.CO ==  0; CONAME = 'LB'; end;
TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,'$\eta_T=',num2str(MODEL.eta_T),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
TITLE = [TITLE, '$(P,J)=(',num2str(GRID.pmaxe),',',num2str(GRID.jmaxe),')$'];