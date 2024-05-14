function [ TOPLOT ] = process_field( DATA, OPTIONS )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Setup time
%%
FRAMES = zeros(size(OPTIONS.TIME));
for i = 1:numel(OPTIONS.TIME)
    [~,FRAMES(i)] =min(abs(OPTIONS.TIME(i)-DATA.Ts3D));
end
FRAMES = unique(FRAMES);
TIME   = DATA.Ts3D(FRAMES);
%% Setup the plot geometry
[KX, KY] = meshgrid(DATA.grids.kx, DATA.grids.ky);
directions = {'y','x','z'};
Nx = DATA.grids.Nx; Ny = DATA.grids.Ny; Nz = DATA.grids.Nz; Nt = numel(FRAMES);
LTXNAME = OPTIONS.NAME;
SKIP_COMP = 0; % Turned on only for kin. distr. func. plot;
Z = 0;
switch OPTIONS.PLAN
    case 'xy'
        XNAME = '$x/\rho_s$'; YNAME = '$y/\rho_s$';
        [X,Y] = meshgrid(DATA.grids.x,DATA.grids.y);
        REALP = 1; COMPDIM = 3; POLARPLOT = 0; SCALE = 1;
    case 'xz'
        XNAME = '$x/\rho_s$'; YNAME = '$z/L_\parallel$';
        [Y,X] = meshgrid(DATA.grids.z,DATA.grids.x);
        REALP = 1; COMPDIM = 1; SCALE = 0; POLARPLOT = 0;
    case 'yz'
        XNAME = '$y/\rho_s$'; YNAME = '$z/L_\parallel$'; 
        [Y,X] = meshgrid(DATA.grids.z,DATA.grids.y);
        REALP = 1; COMPDIM = 2; POLARPLOT = 0; SCALE = 0;
    case 'kxky'
        XNAME = '$k_x\rho_s$'; YNAME = '$k_y\rho_s$';
        [X,Y] = meshgrid(DATA.grids.kx,DATA.grids.ky);
        REALP = 0; COMPDIM = 3; POLARPLOT = 0; SCALE = 1;
    case 'kxz'
        XNAME = '$k_x\rho_s$'; YNAME = '$z/L_\parallel$';
        [Y,X] = meshgrid(DATA.grids.z,DATA.grids.kx);
        REALP = 0; COMPDIM = 1; POLARPLOT = 0; SCALE = 0;
    case 'kyz'
        XNAME = '$k_y\rho_s$'; YNAME = '$z/L_\parallel$';
        [Y,X] = meshgrid(DATA.grids.z,DATA.grids.ky);
        REALP = 0; COMPDIM = 2; POLARPLOT = 0; SCALE = 0;
    case 'xyz'
        XNAME = '$x/\rho_s$'; YNAME = '$y/\rho_s$';
        [X,Y] = meshgrid(DATA.grids.x,DATA.grids.y);
        REALP = 1; COMPDIM = 3; POLARPLOT = 0; SCALE = 1;
        SKIP_COMP = 1; OPTIONS.COMP = 3;
end
Xmax_ = max(max(abs(X))); Ymax_ = max(max(abs(Y)));
Lmin_ = min([Xmax_,Ymax_]);
Rx    = Xmax_/Lmin_ * SCALE + (1-SCALE)*1.2; 
Ry    = Ymax_/Lmin_ * SCALE + (1-SCALE)*1.2; 
ASPECT     = [Rx, Ry, 1];
sz = size(X);
FIELD = zeros(sz(1),sz(2),Nt);
%% Process the field to plot
% --
switch OPTIONS.COMP
    case 'avg'
        compr = @(x) mean(x,COMPDIM);
        COMPNAME = ['avg',directions{COMPDIM}];
        FIELDNAME = ['\langle ',LTXNAME,'\rangle_',directions{COMPDIM}];
    case 'max'
        compr = @(x) max(x,[],COMPDIM);
        COMPNAME = ['max',directions{COMPDIM}];
        FIELDNAME = ['\max_',directions{COMPDIM},' ',LTXNAME];
    otherwise
    switch COMPDIM
        case 1
            i = OPTIONS.COMP;
            compr = @(x) x(i,:,:);
            if REALP
                COMPNAME = sprintf(['y=','%2.1f'],DATA.grids.y(i));
            else
                COMPNAME = sprintf(['k_y=','%2.1f'],DATA.grids.ky(i));
            end
            FIELDNAME = [LTXNAME,'(',COMPNAME,')'];
        case 2
            i = OPTIONS.COMP;
            compr = @(x) x(:,i,:);
            if REALP
                COMPNAME = sprintf(['x=','%2.1f'],DATA.grids.x(i));
            else
                COMPNAME = sprintf(['k_x=','%2.1f'],DATA.grids.kx(i));
            end
            FIELDNAME = [LTXNAME,'(',COMPNAME,')'];
        case 3
            i = OPTIONS.COMP;
            compr = @(x) x(:,:,i);
            COMPNAME = sprintf(['z=','%2.1f'],DATA.grids.z(i));
            FIELDNAME = [LTXNAME,'(',COMPNAME,')'];
    end
end

switch REALP
    case 1 % Real space plot
        INTERP = OPTIONS.INTERP;
        process = @(x) real(fftshift(ifourier_GENE(x)));
        % process = @(x) real(fftshift(ifft2(x,sz(1),sz(2))));
        shift_x = @(x) x;
        shift_y = @(x) x;
    case 0 % Frequencies plot
        INTERP = 0;
        switch COMPDIM
            case 1
                process = @(x) abs(fftshift(x,2));
                shift_x = @(x) fftshift(x,1);
                shift_y = @(x) fftshift(x,1);
            case 2
                process = @(x) abs((x));
                shift_x = @(x) (x);
                shift_y = @(x) (x);        
            case 3
                process = @(x) abs(fftshift(x,2));
                shift_x = @(x) fftshift(x,2);
                shift_y = @(x) fftshift(x,2); 
        end
end

switch OPTIONS.NAME
    case '\phi' %ES pot
        NAME = 'phi';
        FLD_ = DATA.PHI(:,:,:,FRAMES); % data to plot
        OPE_ = 1;        % Operation on data
    case '\psi' %EM pot
        NAME = 'psi';
        FLD_ = DATA.PSI(:,:,:,FRAMES);
        OPE_ = 1;
    case '\phi^{NZ}' % non-zonal ES pot
        NAME = 'phiNZ';
        FLD_ = DATA.PHI(:,:,:,FRAMES);
        OPE_ = (KY~=0);  
   case 'v_{Ey}' % y-comp of ExB velocity
        NAME = 'vy';
        FLD_ = DATA.PHI(:,:,:,FRAMES);
        OPE_ = -1i*KX;  
   case 'v_{Ex}' % x-comp of ExB velocity
        NAME = 'vx';
        FLD_ = DATA.PHI(:,:,:,FRAMES);
        OPE_ = -1i*KY;  
   case 's_{Ey}' % y-comp of ExB shear
        NAME     = 'sy';
        FLD_ = DATA.PHI(:,:,:,FRAMES);
        OPE_ = KX.^2; 
   case 's_{Ex}' % x-comp of ExB shear
        NAME = 'sx';
        FLD_ = DATA.PHI(:,:,:,FRAMES);
        OPE_ = KY.^2; 
   case 'w_{Ez}' % x-comp of ExB shear
        NAME = 'wEz';
        FLD_ = DATA.PHI(:,:,:,FRAMES);
        OPE_ = (KX.^2+KY.^2).*(abs(KX)<=2).*(abs(KY)<=2); 
   case '\omega_z' % ES pot vorticity
        NAME = 'vorticity';
        FLD_ = DATA.PHI(:,:,:,FRAMES);
        OPE_ = -(KX.^2+KY.^2);        
    case 'N_i^{00}' % first ion gyromoment 
        NAME = 'Ni00';
        FLD_ = DATA.Ni00(:,:,:,FRAMES);
        OPE_ = 1;
    case 'N_e^{00}' % first electron gyromoment 
        NAME = 'Ne00';
        FLD_ = DATA.Ne00(:,:,:,FRAMES);
        OPE_ = 1;
    case 'N_i^{00}-N_e^{00}' % first electron gyromoment 
        NAME = 'Ni00-Ne00';
        FLD_ = (DATA.Ni00(:,:,:,FRAMES)-DATA.Ne00(:,:,:,FRAMES))./(poisson_op+(poisson_op==0));
        OPE_ = 1;
    case 'n_i' % ion density
        NAME = 'ni';
        FLD_ = DATA.DENS_I(:,:,:,FRAMES);
        OPE_ = 1;  
    case 'n_e' % electron density
        NAME = 'ne';
        FLD_ = DATA.DENS_E(:,:,:,FRAMES);
        OPE_ = 1;
    case 'k^2n_e' % electron vorticity
        NAME = 'k2ne';
        FLD_ = DATA.DENS_E(:,:,:,FRAMES);
        OPE_ = -(KX.^2+KY.^2);    
    case 'n_i-n_e' % polarisation
        NAME = 'pol';
        OPE_ = 1;
        FLD_ = DATA.DENS_I(:,:,:,FRAMES)...
              -DATA.DENS_E(:,:,:,FRAMES);
    case 'upar_i' % ion density
        NAME = 'upari';
        FLD_ = DATA.UPAR_I(:,:,:,FRAMES);
        OPE_ = 1;  
    case 'upar_e' % electron density
        NAME = 'upare';
        FLD_ = DATA.UPAR_E(:,:,:,FRAMES);
        OPE_ = 1;
    case 'uper_i' % ion density
        NAME = 'uperi';
        FLD_ = DATA.UPER_I(:,:,:,FRAMES);
        OPE_ = 1;  
    case 'uper_e' % electron density
        NAME = 'upere';
        FLD_ = DATA.UPER_E(:,:,:,FRAMES);
        OPE_ = 1;
    case 'T_i' % ion temperature
        NAME = 'Ti';
        FLD_ = DATA.TEMP_I(:,:,:,FRAMES);
        OPE_ = 1; 
    case 'T_e' % ion temperature
        NAME = 'Te';
        FLD_ = DATA.TEMP_E(:,:,:,FRAMES);
        OPE_ = 1; 
    case 'G_x' % ion particle flux
        NAME = 'Gx';
        FLD_ = 0.*DATA.PHI(:,:,:,FRAMES);
        OPE_ = 1;   
        for it = 1:numel(FRAMES)
            tmp = zeros(DATA.grids.Ny,DATA.grids.Nx,Nz);
            for iz = 1:DATA.grids.Nz
                vx_ = real((ifourier_GENE(-1i*KY.*(DATA.PHI   (:,:,iz,FRAMES(it))))));
                ni_ = real((ifourier_GENE(         DATA.DENS_I(:,:,iz,FRAMES(it)))));
                gx_ = vx_.*ni_;
%                 tmp(:,:,iz) = abs(fftshift((squeeze(fft2(gx_,DATA.grids.Ny,DATA.grids.Nx))),2));
                tmp(:,:,iz) = abs((squeeze(fft2(gx_,DATA.grids.Ny,DATA.grids.Nx))));
            end
            FLD_(:,:,:,it)= squeeze((tmp(1:DATA.grids.Nky,1:DATA.grids.Nkx,:)));
        end  
    case 'n_i T_i' % ion heat flux
        NAME = 'niTi';
        FLD_ = 0.*DATA.PHI(:,:,:,FRAMES);
        OPE_ = 1;   
        for it = 1:numel(FRAMES)
            qx_c = zeros(DATA.grids.Ny,DATA.grids.Nx,Nz);
            for iz = 1:DATA.grids.Nz
                ni_  = real((ifourier_GENE(         DATA.DENS_I(:,:,iz,FRAMES(it)))));
                Ti_  = real((ifourier_GENE(         DATA.TEMP_I(:,:,iz,FRAMES(it)))));
                cx_r = ni_.*Ti_;
                cx_ = fft(cx_r,[],1);
                cx_c(:,:,iz) = fft(cx_ ,[],2);
            end
            FLD_(:,:,:,it)= squeeze(cx_c(1:DATA.grids.Nky,1:DATA.grids.Nkx,:))./numel(qx_c);
        end
    case 'Q_{xi}' % ion heat flux
        NAME = 'Qxi';
        FLD_ = 0.*DATA.PHI(:,:,:,FRAMES);
        OPE_ = 1;   
        for it = 1:numel(FRAMES)
            qx_c = zeros(DATA.grids.Ny,DATA.grids.Nx,Nz);
            for iz = 1:DATA.grids.Nz
                vx_  = real((ifourier_GENE(-1i*KY.*(DATA.PHI   (:,:,iz,FRAMES(it))))));
                ni_  = real((ifourier_GENE(         DATA.DENS_I(:,:,iz,FRAMES(it)))));
                Ti_  = real((ifourier_GENE(         DATA.TEMP_I(:,:,iz,FRAMES(it)))));
                qx_r = vx_.*ni_.*Ti_;
                qx_ = fft(qx_r,[],1);
                qx_c(:,:,iz) = fft(qx_ ,[],2);
            end
            FLD_(:,:,:,it)= squeeze(qx_c(1:DATA.grids.Nky,1:DATA.grids.Nkx,:))./numel(qx_c);
        end
    case 'Q_{xe}' % electron heat flux
        NAME = 'Qxe';
        FLD_ = 0.*DATA.PHI(:,:,:,FRAMES);
        OPE_ = 1;   
        for it = 1:numel(FRAMES)
            qx_c = zeros(DATA.grids.Ny,DATA.grids.Nx,Nz);
            for iz = 1:DATA.grids.Nz
                vx_  = real((ifourier_GENE(-1i*KY.*(DATA.PHI   (:,:,iz,FRAMES(it))))));
                ne_  = real((ifourier_GENE(         DATA.DENS_E(:,:,iz,FRAMES(it)))));
                Te_  = real((ifourier_GENE(         DATA.TEMP_E(:,:,iz,FRAMES(it)))));
                qx_r = vx_.*ne_.*Te_;
                qx_ = fft(qx_r,[],1);
                qx_c(:,:,iz) = fft(qx_ ,[],2);
            end
            FLD_(:,:,:,it)= squeeze(qx_c(1:DATA.grids.Nky,1:DATA.grids.Nkx,:))./numel(qx_c);
        end    
    case 'f_i'
        SKIP_COMP = 1;
        shift_x = @(x) x; shift_y = @(y) y;
        NAME = 'fi'; OPTIONS.SPECIE = 'i';
        for it = 1:numel(FRAMES)
            OPTIONS.T = DATA.Ts5D(FRAMES(it));
            OPTIONS.Z = OPTIONS.COMP;
            [~,~,FIELD(:,:,it)] = compute_fa(DATA,OPTIONS);
        end  
    case 'f_e'
        SKIP_COMP = 1;
        shift_x = @(x) x; shift_y = @(y) y;
        NAME = 'fe'; OPTIONS.SPECIE = 'e';
        [~,it0_] =min(abs(OPTIONS.TIME(1)-DATA.Ts5D));
        [~,it1_]=min(abs(OPTIONS.TIME(end)-DATA.Ts5D));
        dit_ = 1;%ceil((it1_-it0_+1)/10); 
        FRAMES = it0_:dit_:it1_;
        iz = 1;
        for it = 1:numel(FRAMES)
            OPTIONS.T = DATA.Ts5D(FRAMES(it));
            [~,~,FIELD(:,:,it)] = compute_fa(DATA,OPTIONS);
        end  
    otherwise
        disp('Fieldname not recognized');
        return
end
% Process the field according to the 2D plane and the space (real/cpx)
if ~SKIP_COMP
    if COMPDIM == 3
        for it = 1:numel(FRAMES)
            tmp = squeeze(compr(OPE_.*FLD_(:,:,:,it)));
            FIELD(:,:,it) = squeeze(process(tmp));
        end
    else
        if REALP
            tmp = zeros(Ny,Nx,Nz);
        else
            tmp = zeros(DATA.grids.Nky,DATA.grids.Nkx,Nz);
        end
        for it = 1:numel(FRAMES)
            for iz = 1:numel(DATA.grids.z)
                tmp(:,:,iz) = squeeze(process(OPE_.*FLD_(:,:,iz,it)));
            end
            FIELD(:,:,it) = squeeze(compr(tmp));
        end                
    end
else
    clear FIELD;
    FIELD = squeeze(process(FLD_));
end
TOPLOT.FIELD     = FIELD;
TOPLOT.TIME      = TIME;
TOPLOT.FRAMES    = FRAMES;
TOPLOT.X         = shift_x(X);
TOPLOT.Y         = shift_y(Y);
TOPLOT.Z         = Z;
TOPLOT.FIELDNAME = FIELDNAME;
TOPLOT.XNAME     = XNAME;
TOPLOT.YNAME     = YNAME;
TOPLOT.FILENAME  = [NAME,'_',OPTIONS.PLAN,'_',COMPNAME];
TOPLOT.ASPECT    = ASPECT;
TOPLOT.FRAMES    = FRAMES;
TOPLOT.INTERP    = INTERP;

end

