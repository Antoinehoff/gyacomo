addpath(genpath('../matlab')) % ... add
%% Grid configuration
N       = 200;     % Frequency gridpoints (Nkr = N/2)
L       = 120;     % Size of the squared frequency domain
dk      = 2.*pi/L;
kmax    = N/2*dk;
kr      = dk*(0:N/2);
kz      = dk*(0:N/2);
[KZ, KR]= meshgrid(kz,kr);
KPERP   = sqrt(KR.^2 + KZ.^2);
kperp   = reshape(KPERP,[1,numel(kr)^2]);
kperp   = uniquetol(kperp,1e-14);
Nperp   = numel(kperp);
%% Model
    % ! 0 = nothing 
    % ! 1 = coulomb
    % ! 2 = pitch-angle (with constant coll.freq.)
    % ! 3 = sugama
    % ! 4 = pitch-angle with veloty dependent coll. freq.
    % ! 5 = improved sugama
    % ! 6 = Hirschman-Sigmar Clarke
    % ! 7 = Abel (only for like species)
    % ! 8 = conserving momentun pitch-angle operator (with veloty dependent coll. freq.)
    % ! 9 = GK coulomb polarization term
CO = 3;
GK = 1;
P = 10;
J = 5;
M_FLR= 0; %to increase the NFLR
if 0
    %% plot the kperp distribution
   figure
   plot(kperp)
end
% %% Check if the differences btw kperp is larger than naming precision
%%
%% We compute only on a kperp grid with dk space from 0 to kperpmax
kperpmax = sqrt(2) * kmax;
kperp = unique([0:dk:kperpmax,kperpmax]);
%% Naming
if CO == 1; CONAME = 'FC'; end;
if CO == 2; CONAME = 'PA'; end;
if CO == 3; CONAME = 'SG'; end;
matfilename = ['../iCa/',CONAME,'_P_',num2str(P),'_J_',num2str(J),...
    '_N_',num2str(N),'_dk_',num2str(dk),'_MFLR_',num2str(M_FLR),'.h5'];

n_ = 1;
for k_ = kperp
    disp(['-Writing matrix for kperp = ',num2str(k_)])
    %% Script to run COSOlver in order to create needed collision matrices
    COSOlver_path = '../../Documents/MoliSolver/COSOlver/';
    COSOLVER.pmaxe = P;
    COSOLVER.jmaxe = J;
    COSOLVER.pmaxi = P;
    COSOLVER.jmaxi = J;
    COSOLVER.kperp = k_;

    COSOLVER.neFLR    = min(ceil((2/3*kperpmax)^2),max(5,ceil(COSOLVER.kperp^2)))+M_FLR; % rule of thumb for sum truncation
    COSOLVER.niFLR    = min(ceil((2/3*kperpmax)^2),max(5,ceil(COSOLVER.kperp^2)))+M_FLR;
    COSOLVER.idxT4max = 40;

    COSOLVER.neFLRs = 0; %  ... only for GK abel 
    COSOLVER.npeFLR = 0; %  ... only for GK abel 
    COSOLVER.niFLRs = 0; %  ... only for GK abel 
    COSOLVER.npiFLR = 0; %  ... only for GK abel 

    COSOLVER.gk = GK; 
    COSOLVER.co = CO;
    
    k_string      = sprintf('%0.4f',k_);
    self_mat_file_name = ['../iCa/self_Coll_GKE_',num2str(COSOLVER.gk),'_GKI_',num2str(COSOLVER.gk),...
        '_ESELF_',num2str(COSOLVER.co),'_ISELF_',num2str(COSOLVER.co),...
        '_Pmaxe_',num2str(COSOLVER.pmaxe),'_Jmaxe_',num2str(COSOLVER.jmaxe),...
        '_Pmaxi_',num2str(COSOLVER.pmaxi),'_Jmaxi_',num2str(COSOLVER.jmaxi),...
        '_JE_12',...
        '_NFLR_',num2str(COSOLVER.neFLR),...
        '_kperp_',k_string,'.h5'];
    ei_mat_file_name = ['../iCa/ei_Coll_GKE_',num2str(COSOLVER.gk),'_GKI_',num2str(COSOLVER.gk),...
        '_ETEST_',num2str(COSOLVER.co),'_EBACK_',num2str(COSOLVER.co),...
        '_Pmaxe_',num2str(COSOLVER.pmaxe),'_Jmaxe_',num2str(COSOLVER.jmaxe),...
        '_Pmaxi_',num2str(COSOLVER.pmaxi),'_Jmaxi_',num2str(COSOLVER.jmaxi),...
        '_JE_12_tau_1.0000_mu_0.0233',...
        '_NFLRe_',num2str(COSOLVER.neFLR),'_NFLRi_',num2str(COSOLVER.neFLR)...
        '_kperp_',k_string,'.h5'];
    ie_mat_file_name = ['../iCa/ie_Coll_GKE_',num2str(COSOLVER.gk),'_GKI_',num2str(COSOLVER.gk),...
        '_ITEST_',num2str(COSOLVER.co),'_IBACK_',num2str(COSOLVER.co),...
        '_Pmaxe_',num2str(COSOLVER.pmaxe),'_Jmaxe_',num2str(COSOLVER.jmaxe),...
        '_Pmaxi_',num2str(COSOLVER.pmaxi),'_Jmaxi_',num2str(COSOLVER.jmaxi),...
        '_JE_12_tau_1.0000_mu_0.0233',...
        '_NFLRe_',num2str(COSOLVER.neFLR),'_NFLRi_',num2str(COSOLVER.neFLR)...
        '_kperp_',k_string,'.h5'];
   
    %% Load self matrix
    C_self = h5read(self_mat_file_name,'/Caapj/Ceepj');
    sz_ = size(C_self);
    % Write it in the compiled file
    h5create(matfilename,['/Caapj/',k_string],sz_)
    h5write(matfilename,['/Caapj/',k_string],C_self)
    
    %% Load ei matrices
    % Field
    C_eiF = h5read(ei_mat_file_name,'/Ceipj/CeipjF');
    sz_ = size(C_eiF);
    h5create(matfilename,['/CeipjF/',k_string],sz_)
    h5write(matfilename,['/CeipjF/',k_string],C_eiF)
    % Test
    C_eiT = h5read(ei_mat_file_name,'/Ceipj/CeipjT');
    sz_ = size(C_eiT);
    h5create(matfilename,['/CeipjT/',k_string],sz_)
    h5write(matfilename,['/CeipjT/',k_string],C_eiT)
    
    %% Load ie matrices
    % Field
    C_ieF = h5read(ie_mat_file_name,'/Ciepj/CiepjF');
    sz_ = size(C_ieF);
    h5create(matfilename,['/CiepjF/',k_string],sz_)
    h5write(matfilename,['/CiepjF/',k_string],C_ieF)
    % Test
    C_ieT = h5read(ie_mat_file_name,'/Ciepj/CiepjT');
    sz_ = size(C_eiT);
    h5create(matfilename,['/CiepjT/',k_string],sz_)
    h5write(matfilename,['/CiepjT/',k_string],C_ieT)
    
    %% Copy fort.90 input file and put grid params
    if(k_ == 0)
        h5create(matfilename,'/dk',1);
        h5write (matfilename,'/dk',dk);   
        h5create(matfilename,'/N',1);
        h5write (matfilename,'/N',N);
        h5create(matfilename,'/Pmaxe',1);
        h5write (matfilename,'/Pmaxe',P);   
        h5create(matfilename,'/Jmaxe',1);
        h5write (matfilename,'/Jmaxe',J);   
        h5create(matfilename,'/Pmaxi',1);
        h5write (matfilename,'/Pmaxi',P);  
        h5create(matfilename,'/Jmaxi',1);
        h5write (matfilename,'/Jmaxi',J);   
    end
    
end
disp(['File saved @ :',matfilename])