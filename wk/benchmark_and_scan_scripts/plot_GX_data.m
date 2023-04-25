figure
% Linear CBC adiab. e, miller (figure 1 top)
GX_data_ = load('/home/ahoffman/GX_paper_data/fig_1_1.txt');
subplot(2,2,1)
scale_x = 1.0; 
scale_y = 2.777778; % gamma_AH/gamma_NM = R/a
plot(scale_x*GX_data_(:,1),scale_y*GX_data_(:,2),'--ob','Displayname','GX');
xlabel('$k_y\rho_i$'); ylabel('$\gamma R/v_{ti}$');
title('(figure 1 top)');
% Linear CBC kin. e, miller (figure 2 top left)
GX_data_ = load('/home/ahoffman/GX_paper_data/fig_2_1.txt');
subplot(2,2,2)
scale_x = 1.0; 
scale_y = 2.777778; % gamma_AH/gamma_NM = R/a
plot(scale_x*GX_data_(:,1),scale_y*GX_data_(:,2),'--ob','Displayname','GX');
xlabel('$k_y\rho_i$'); ylabel('$\gamma R/v_{ti}$');
title('(figure 2 top left)');
% Nonlinear heatflux, adiab. e, miller (figure 4 left)
GX_data_ = load('/home/ahoffman/GX_paper_data/fig_4_1.txt');
subplot(2,2,3)
scale_x = 1/2.777778;  % t_AH/t_NM = (a/R)
scale_y = 2.777778^2;
plot(scale_x*GX_data_(:,1),scale_y*GX_data_(:,2),'-b','Displayname','GX');
xlabel('$tv_{ti}/R$'); ylabel('$Q_i/G_{GB}^R$');
title('(figure 4 right)');
% Dimits shift, adiab. e, miller (figure 4 right)
GX_data_ = load('/home/ahoffman/GX_paper_data/fig_4_2.txt');
subplot(2,2,4)
scale_x = 2.777778; % wT_AH/wT_NM = a/R
scale_y = 2.777778^2; % Q_AH/Q_NM = (a/R)^2
plot(scale_x*GX_data_(:,1),scale_y*GX_data_(:,2),'--ob','Displayname','GX');
xlabel('$R/L_T$'); ylabel('$Q_i/Q_{GB}^R$');
title('(figure 4 right)');