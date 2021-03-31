addpath(genpath('../matlab')) % ... add
%% Grid configuration
N       = 10;     % Frequency gridpoints (Nkr = N/2)
L       = 120;     % Size of the squared frequency domain
dk      = 2*pi/L;
kmax    = N/2*dk;
kr      = dk*(0:N/2);
kz      = dk*(0:N/2);
[KZ, KR]= meshgrid(kz,kr);
KPERP   = sqrt(KR.^2 + KZ.^2);
kperp   = reshape(KPERP,[1,numel(kr)^2]);
kperp   = uniquetol(kperp,1e-14);
Nperp   = numel(kperp);
COSOlver_path = '../../Documents/MoliSolver/COSOlver/';
COSOLVER.pmaxe = 10;
COSOLVER.jmaxe = 5;
COSOLVER.pmaxi = 10;
COSOLVER.jmaxi = 5;

COSOLVER.neFLR  = max(5,ceil(COSOLVER.kperp^2)); % rule of thumb for sum truncation
COSOLVER.niFLR  = max(5,ceil(COSOLVER.kperp^2));

COSOLVER.neFLRs = 0; %  ... only for GK abel 
COSOLVER.npeFLR = 0; %  ... only for GK abel 
COSOLVER.niFLRs = 0; %  ... only for GK abel 
COSOLVER.npiFLR = 0; %  ... only for GK abel 

COSOLVER.gk = 1; 
COSOLVER.co = 3;
if 0
    %% plot the kperp distribution
   figure
   plot(kperp)
end
%% Load the matrices
C_self_i = zeros((COSOLVER.pmaxi+1)*(COSOLVER.jmaxi+1),(COSOLVER.pmaxi+1)*(COSOLVER.jmaxi+1),Nperp);

for n_ = 1:Nperp
    COSOLVER.kperp = kperp(n_);
    k_string      = sprintf('%0.4f',kperp(n_));
    mat_file_name = ['../iCa/self_Coll_GKE_',num2str(COSOLVER.gk),'_GKI_',num2str(COSOLVER.gk),...
    '_ESELF_',num2str(COSOLVER.co),'_ISELF_',num2str(COSOLVER.co),...
    '_Pmaxe_',num2str(COSOLVER.pmaxe),'_Jmaxe_',num2str(COSOLVER.jmaxe),...
    '_Pmaxi_',num2str(COSOLVER.pmaxi),'_Jmaxi_',num2str(COSOLVER.jmaxi),...
    '_JE_12',...
    '_NFLR_',num2str(COSOLVER.neFLR),...
    '_kperp_',k_string,'.h5'];

    tmp  = h5read(mat_file_name,'/Caapj/Ceepj');
    C_self_i(:,:,n_) = tmp;
end

%% Post processing
dC_self_i = diff(C_self_i,1,3);
gvar_dC    = squeeze(sum(sum(abs(dC_self_i),2),1));
%% Plots
%% all coeffs evolution
figure
for ip_ = 1:COSOLVER.pmaxi+1
    for ij_ = 1:COSOLVER.jmaxi+1
        plot(kperp,squeeze(C_self_i(ip_,ij_,:)),'o'); hold on;
    end
end
%% global matrix variation
figure;
kperpmid = 0.5*(kperp(2:end)+kperp(1:end-1));
plot(kperpmid,gvar_dC./diff(kperp)');





