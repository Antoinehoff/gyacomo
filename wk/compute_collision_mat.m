addpath(genpath('../matlab')) % ... add
%% Grid configuration
N       = 200;     % Frequency gridpoints (Nkr = N/2)
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
if 0
    %% plot the kperp distribution
   figure
   plot(kperp)
end
% %% Check if the differences btw kperp is larger than naming precision
% dkperp  = diff(kperp);
% warning = sum(dkperp<0.0002);
% if warning > 0
%     disp('Warning : dkperp < 0.0002');
% end
%%
%% We compute only on a kperp grid with dk space from 0 to kperpmax
kperp = unique([0:dk:(sqrt(2)*kmax),sqrt(2)*kmax]);
kperpmax = sqrt(2) * kmax;
%%
n_ = 1;
for k_ = kperp
    disp(['----------Computing matrix ',num2str(n_),'/',num2str(Nperp),'----------'])
    %% Script to run COSOlver in order to create needed collision matrices
    COSOlver_path = '../../Documents/MoliSolver/COSOlver/';
    COSOLVER.pmaxe = 20;
    COSOLVER.jmaxe = 10;
    COSOLVER.pmaxi = 20;
    COSOLVER.jmaxi = 10;
    COSOLVER.kperp = k_;

    COSOLVER.neFLR    = min(ceil((2/3*kperpmax)^2),max(5,ceil(COSOLVER.kperp^2))); % rule of thumb for sum truncation
    COSOLVER.niFLR    = min(ceil((2/3*kperpmax)^2),max(5,ceil(COSOLVER.kperp^2)));
    COSOLVER.idxT4max = 40;

    COSOLVER.neFLRs = 0; %  ... only for GK abel 
    COSOLVER.npeFLR = 0; %  ... only for GK abel 
    COSOLVER.niFLRs = 0; %  ... only for GK abel 
    COSOLVER.npiFLR = 0; %  ... only for GK abel 

    COSOLVER.gk = 1; 
    COSOLVER.co = 3;
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

    write_fort90_COSOlver
    
    k_string      = sprintf('%0.4f',k_);
    mat_file_name = ['../iCa/self_Coll_GKE_',num2str(COSOLVER.gk),'_GKI_',num2str(COSOLVER.gk),...
        '_ESELF_',num2str(COSOLVER.co),'_ISELF_',num2str(COSOLVER.co),...
        '_Pmaxe_',num2str(COSOLVER.pmaxe),'_Jmaxe_',num2str(COSOLVER.jmaxe),...
        '_Pmaxi_',num2str(COSOLVER.pmaxi),'_Jmaxi_',num2str(COSOLVER.jmaxi),...
        '_JE_12',...
        '_NFLR_',num2str(COSOLVER.neFLR),...
        '_kperp_',k_string,'.h5'];
    if (exist(mat_file_name,'file') > 0)
        disp(['Matrix available for kperp = ',k_string]);
    else
        cd ../../Documents/MoliSolver/COSOlver/
        disp(['Matrix not found for kperp = ',k_string]);
        disp([num2str(n_),'/',Nperp])
        disp('computing...');
        CMD = 'mpirun -np 6 bin/CO 2 2 2 > out.txt';
        disp(CMD); 
        system(CMD);
        system(CMD);
        disp('..done');
        cd ../../../HeLaZ/wk
    end
    n_ = n_ + 1;
end