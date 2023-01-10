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
%% Setup the plot geometry
[KX, KY] = meshgrid(DATA.kx, DATA.ky);
directions = {'x','y','z'};
Nx = DATA.Nx; Ny = DATA.Ny; Nz = DATA.Nz; Nt = numel(FRAMES);
POLARPLOT = OPTIONS.POLARPLOT;
LTXNAME = OPTIONS.NAME;
switch OPTIONS.PLAN
    case 'xy'
        XNAME = '$x$'; YNAME = '$y$';
        [X,Y] = meshgrid(DATA.x,DATA.y);
        REALP = 1; COMPDIM = 3; POLARPLOT = 0; SCALE = 1;
    case 'xz'
        XNAME = '$x$'; YNAME = '$z$';
        [Y,X] = meshgrid(DATA.z,DATA.x);
        REALP = 1; COMPDIM = 1; SCALE = 0;
    case 'yz'
        XNAME = '$y$'; YNAME = '$z$'; 
        [Y,X] = meshgrid(DATA.z,DATA.y);
        REALP = 1; COMPDIM = 2; SCALE = 0;
    case 'kxky'
        XNAME = '$k_x$'; YNAME = '$k_y$';
        [X,Y] = meshgrid(DATA.kx,DATA.ky);
        REALP = 0; COMPDIM = 3; POLARPLOT = 0; SCALE = 1;
    case 'kxz'
        XNAME = '$k_x$'; YNAME = '$z$';
        [Y,X] = meshgrid(DATA.z,DATA.kx);
        REALP = 0; COMPDIM = 1; POLARPLOT = 0; SCALE = 0;
    case 'kyz'
        XNAME = '$k_y$'; YNAME = '$z$';
        [Y,X] = meshgrid(DATA.z,DATA.ky);
        REALP = 0; COMPDIM = 2; POLARPLOT = 0; SCALE = 0;
end
Xmax_ = max(max(abs(X))); Ymax_ = max(max(abs(Y)));
Lmin_ = min([Xmax_,Ymax_]);
Rx    = Xmax_/Lmin_ * SCALE + (1-SCALE)*1.2; 
Ry    = Ymax_/Lmin_ * SCALE + (1-SCALE)*1.2; 
ASPECT     = [Rx, Ry, 1];
DIMENSIONS = [500, 1000, OPTIONS.RESOLUTION*Rx, OPTIONS.RESOLUTION*Ry];
% Polar grid
POLARNAME = [];
if POLARPLOT
    POLARNAME = 'polar';
    X__ = (X+DATA.a).*cos(Y);
    Y__ = (X+DATA.a).*sin(Y);
    X = X__;
    Y = Y__;
    XNAME='X';
    YNAME='Z';
    DIMENSIONS = [100, 100, OPTIONS.RESOLUTION, OPTIONS.RESOLUTION];
    ASPECT     = [1,1,1];
    sz = size(X);
    FIELD = zeros(sz(1),sz(2),Nt);
else
    sz = size(X);
    FIELD = zeros(sz(1),sz(2),Nt);
end
%% Process the field to plot

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
                COMPNAME = sprintf(['y=','%2.1f'],DATA.x(i));
            else
                COMPNAME = sprintf(['k_y=','%2.1f'],DATA.kx(i));
            end
            FIELDNAME = [LTXNAME,'(',COMPNAME,')'];
        case 2
            i = OPTIONS.COMP;
            compr = @(x) x(:,i,:);
            if REALP
                COMPNAME = sprintf(['x=','%2.1f'],DATA.y(i));
            else
                COMPNAME = sprintf(['k_x=','%2.1f'],DATA.ky(i));
            end
            FIELDNAME = [LTXNAME,'(',COMPNAME,')'];
        case 3
            i = OPTIONS.COMP;
            compr = @(x) x(:,:,i);
            COMPNAME = sprintf(['z=','%2.1f'],DATA.z(i));
            FIELDNAME = [LTXNAME,'(',COMPNAME,')'];
    end
end

switch REALP
    case 1 % Real space plot
        INTERP = OPTIONS.INTERP;
        process = @(x) real(fftshift(ifourier_GENE(x)));
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
        FLD_ = DATA.DENS_I(:,:,:,FRAMES)+DATA.DENS_E(:,:,:,FRAMES);
    case 'T_i' % ion temperature
        NAME = 'Ti';
        FLD_ = DATA.TEMP_I(:,:,:,FRAMES);
        OPE_ = 1; 
    case 'G_x' % ion particle flux
        NAME = 'Gx';
        FLD_ = 0.*DATA.PHI(:,:,:,FRAMES);
        OPE_ = 1;   
        for it = 1:numel(FRAMES)
            tmp = zeros(DATA.Ny,DATA.Nx,Nz);
            for iz = 1:DATA.Nz
                vx_ = real((ifourier_GENE(-1i*KY.*(DATA.PHI   (:,:,iz,FRAMES(it))))));
                ni_ = real((ifourier_GENE(         DATA.DENS_I(:,:,iz,FRAMES(it)))));
                gx_ = vx_.*ni_;
%                 tmp(:,:,iz) = abs(fftshift((squeeze(fft2(gx_,DATA.Ny,DATA.Nx))),2));
                tmp(:,:,iz) = abs((squeeze(fft2(gx_,DATA.Ny,DATA.Nx))));
            end
            FLD_(:,:,it)= squeeze(compr(tmp(1:DATA.Nky,1:DATA.Nkx,:)));
        end   
    case 'Q_x' % ion heat flux
        NAME = 'Qx';
        FLD_ = 0.*DATA.PHI(:,:,:,FRAMES);
        OPE_ = 1;   
        for it = 1:numel(FRAMES)
            tmp = zeros(DATA.Ny,DATA.Nx,Nz);
            for iz = 1:DATA.Nz
                vx_ = real((ifourier_GENE(-1i*KY.*(DATA.PHI   (:,:,iz,FRAMES(it))))));
                ni_ = real((ifourier_GENE(         DATA.DENS_I(:,:,iz,FRAMES(it)))));
                Ti_ = real((ifourier_GENE(         DATA.TEMP_I(:,:,iz,FRAMES(it)))));
                qx_ = vx_.*ni_.*Ti_;
%                 tmp(:,:,iz) = abs(fftshift((squeeze(fft2(gx_,DATA.Ny,DATA.Nx))),2));
                tmp(:,:,iz) = abs((squeeze(fft2(qx_,DATA.Ny,DATA.Nx))));
            end
            FLD_(:,:,it)= squeeze(compr(tmp(1:DATA.Nky,1:DATA.Nkx,:)));
        end     
    otherwise
        disp('Fieldname not recognized');
        return
end
% Process the field according to the 2D plane and the space (real/cpx)
if COMPDIM == 3
    for it = 1:numel(FRAMES)
        tmp = squeeze(compr(OPE_.*FLD_(:,:,:,it)));
        FIELD(:,:,it) = squeeze(process(tmp));
    end
else
    if REALP
        tmp = zeros(Ny,Nx,Nz);
    else
        tmp = zeros(DATA.Nky,DATA.Nkx,Nz);
    end
    for it = 1:numel(FRAMES)
        for iz = 1:numel(DATA.z)
            tmp(:,:,iz) = squeeze(process(OPE_.*FLD_(:,:,:,it)));
        end
        FIELD(:,:,it) = squeeze(compr(tmp));
    end                
end
TOPLOT.FIELD     = FIELD;
TOPLOT.FRAMES    = FRAMES;
TOPLOT.X         = shift_x(X);
TOPLOT.Y         = shift_y(Y);
TOPLOT.FIELDNAME = FIELDNAME;
TOPLOT.XNAME     = XNAME;
TOPLOT.YNAME     = YNAME;
TOPLOT.FILENAME  = [NAME,'_',OPTIONS.PLAN,'_',COMPNAME,'_',POLARNAME];
TOPLOT.DIMENSIONS= DIMENSIONS;
TOPLOT.ASPECT    = ASPECT;
TOPLOT.FRAMES    = FRAMES;
TOPLOT.INTERP    = INTERP;
end

