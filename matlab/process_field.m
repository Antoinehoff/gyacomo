function [ TOPLOT ] = process_field( DATA, OPTIONS )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Setup time
%%
FRAMES = zeros(size(OPTIONS.TIME));
for i = 1:numel(OPTIONS.TIME)
    [~,FRAMES(i)] =min(abs(OPTIONS.TIME(i)-DATA.Ts3D));
end
%% Setup the plot geometry
[KY,~] = meshgrid(DATA.ky,DATA.kx);
directions = {'x','y','z'};
Nx = DATA.Nx; Ny = DATA.Ny; Nz = DATA.Nz; Nt = numel(OPTIONS.TIME);
POLARPLOT = OPTIONS.POLARPLOT;
LTXNAME = OPTIONS.NAME;
switch OPTIONS.PLAN
    case 'xy'
        XNAME = '$x$'; YNAME = '$y$';
        [Y,X] = meshgrid(DATA.y,DATA.x);
        REALP = 1; COMPDIM = 3; POLARPLOT = 0; SCALE = 1;
    case 'xz'
        XNAME = '$x$'; YNAME = '$z$';
        [Y,X] = meshgrid(DATA.z,DATA.x);
        REALP = 1; COMPDIM = 2; SCALE = 0;
    case 'yz'
        XNAME = '$y$'; YNAME = '$z$'; 
        [Y,X] = meshgrid(DATA.z,DATA.y);
        REALP = 1; COMPDIM = 1; SCALE = 0;
    case 'kxky'
        XNAME = '$k_x$'; YNAME = '$k_y$';
        [Y,X] = meshgrid(DATA.ky,DATA.kx);
        REALP = 0; COMPDIM = 3; POLARPLOT = 0; SCALE = 1;
    case 'kxz'
        XNAME = '$k_x$'; YNAME = '$z$';
        [Y,X] = meshgrid(DATA.z,DATA.kx);
        REALP = 0; COMPDIM = 2; POLARPLOT = 0; SCALE = 0;
    case 'kyz'
        XNAME = '$k_y$'; YNAME = '$z$';
        [Y,X] = meshgrid(DATA.z,DATA.ky);
        REALP = 0; COMPDIM = 1; POLARPLOT = 0; SCALE = 0;
end
Xmax_ = max(max(abs(X))); Ymax_ = max(max(abs(Y)));
Lmin_ = min([Xmax_,Ymax_]);
Rx    = Xmax_/Lmin_ * SCALE + (1-SCALE)*1.2; 
Ry    = Ymax_/Lmin_ * SCALE + (1-SCALE)*1.2; 
DIMENSIONS = [100, 100, 400*Rx, 400*Ry];
ASPECT     = [Rx, Ry, 1];
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
    DIMENSIONS = [100, 100, 500, 500];
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
                COMPNAME = sprintf(['x=','%2.1f'],DATA.x(i));
            else
                COMPNAME = sprintf(['k_x=','%2.1f'],DATA.kx(i));
            end
            FIELDNAME = [LTXNAME,'(',COMPNAME,')'];
        case 2
            i = OPTIONS.COMP;
            compr = @(x) x(:,i,:);
            if REALP
                COMPNAME = sprintf(['y=','%2.1f'],DATA.y(i));
            else
                COMPNAME = sprintf(['k_y=','%2.1f'],DATA.ky(i));
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
        process = @(x) real(fftshift(ifft2(x,Nx,Ny)));
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
    case '\phi'
        NAME = 'phi';
        if COMPDIM == 3
            for it = FRAMES
                tmp = squeeze(compr(DATA.PHI(:,:,:,it)));
                FIELD(:,:,it) = squeeze(process(tmp));
            end
        else
            if REALP
                tmp = zeros(Nx,Ny,Nz);
            else
                tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            end
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(DATA.PHI(:,:,iz,it)));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end                
        end
    case 'n_e'
        NAME = 'ne';
        if COMPDIM == 3
            for it = FRAMES
                tmp = squeeze(compr(DATA.DENS_E(:,:,:,it)));
                FIELD(:,:,it) = squeeze(process(tmp));
            end
        else
            if REALP
                tmp = zeros(Nx,Ny,Nz);
            else
                tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            end
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(DATA.DENS_E(:,:,iz,it)));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end                
        end
    case 'k^2n_e'
        NAME = 'k2ne';
        [KY, KX] = meshgrid(DATA.ky, DATA.kx);
        if COMPDIM == 3
            for it = FRAMES
                for iz = 1:DATA.Nz
                tmp = squeeze(compr(-(KX.^2+KY.^2).*DATA.DENS_E(:,:,iz,it)));
                end
                FIELD(:,:,it) = squeeze(process(tmp));
            end
        else
            if REALP
                tmp = zeros(Nx,Ny,Nz);
            else
                tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            end
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(-(KX.^2+KY.^2).*DATA.DENS_E(:,:,iz,it)));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end                
        end  
    case 'n_e^{NZ}'
        NAME = 'neNZ';
        if COMPDIM == 3
            for it = FRAMES
                tmp = squeeze(compr(DATA.DENS_E(:,:,:,it).*(KY~=0)));
                FIELD(:,:,it) = squeeze(process(tmp));
            end
        else
            if REALP
                tmp = zeros(Nx,Ny,Nz);
            else
                tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            end
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(DATA.DENS_E(:,:,iz,it).*(KY~=0)));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end                
        end        
    case 'n_i'
        NAME = 'ni';
        if COMPDIM == 3
            for it = FRAMES
                tmp = squeeze(compr(DATA.DENS_I(:,:,:,it)));
                FIELD(:,:,it) = squeeze(process(tmp));
            end
        else
            if REALP
                tmp = zeros(Nx,Ny,Nz);
            else
                tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            end
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(DATA.DENS_I(:,:,iz,it)));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end                
        end
    case 'n_i^{NZ}'
        NAME = 'niNZ';
        if COMPDIM == 3
            for it = FRAMES
                tmp = squeeze(compr(DATA.DENS_I(:,:,:,it).*(KY~=0)));
                FIELD(:,:,it) = squeeze(process(tmp));
            end
        else
            if REALP
                tmp = zeros(Nx,Ny,Nz);
            else
                tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            end
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(DATA.DENS_I(:,:,iz,it).*(KY~=0)));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end                
        end
    case '\phi^{NZ}'
        NAME = 'phiNZ';
        if COMPDIM == 3
            for it = FRAMES
                tmp = squeeze(compr(DATA.PHI(:,:,:,it).*(KY~=0)));
                FIELD(:,:,it) = squeeze(process(tmp));
            end
        else
            if REALP
                tmp = zeros(Nx,Ny,Nz);
            else
                tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            end
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(DATA.PHI(:,:,iz,it).*(KY~=0)));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end                
        end
   case 'v_y'
        NAME     = 'vy';
        [~, KX] = meshgrid(DATA.ky, DATA.kx);
        vE      = zeros(DATA.Nx,DATA.Ny,DATA.Nz,numel(FRAMES));
        for it = FRAMES % Compute the 3D real transport for each frame
            for iz = 1:DATA.Nz
            vE(:,:,iz,it)  = real(ifft2(-1i*KX.*(DATA.PHI(:,:,iz,it)),DATA.Nx,DATA.Ny));
            end
        end
        if REALP % plot in real space
            for it = FRAMES
                FIELD(:,:,it) = squeeze(compr(vE(:,:,:,it)));
            end
        else % Plot spectrum
            process = @(x) abs(fftshift(x,2));
            tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(-1i*KX.*(DATA.PHI(:,:,iz,it))));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end   
        end 
   case 'v_x'
        NAME     = 'vx';
        [KY, ~] = meshgrid(DATA.ky, DATA.kx);
        vE      = zeros(DATA.Nx,DATA.Ny,DATA.Nz,numel(FRAMES));
        for it = FRAMES % Compute the 3D real transport for each frame
            for iz = 1:DATA.Nz
            vE(:,:,iz,it)  = real(ifft2(-1i*KY.*(DATA.PHI(:,:,iz,it)),DATA.Nx,DATA.Ny));
            end
        end
        if REALP % plot in real space
            for it = FRAMES
                FIELD(:,:,it) = squeeze(compr(vE(:,:,:,it)));
            end
        else % Plot spectrum
            process = @(x) abs(fftshift(x,2));
            tmp = zeros(DATA.Nkx,DATA.Nky,Nz);
            for it = FRAMES
                for iz = 1:numel(DATA.z)
                    tmp(:,:,iz) = squeeze(process(-1i*KY.*(DATA.PHI(:,:,iz,it))));
                end
                FIELD(:,:,it) = squeeze(compr(tmp));
            end   
        end 
    case '\Gamma_x'
    NAME     = 'Gx';
    [KY, ~] = meshgrid(DATA.ky, DATA.kx);
    Gx      = zeros(DATA.Nx,DATA.Ny,DATA.Nz,numel(FRAMES));
    for it = FRAMES % Compute the 3D real transport for each frame
        for iz = 1:DATA.Nz
        Gx(:,:,iz,it)  = real((ifft2(-1i*KY.*(DATA.PHI(:,:,iz,it)),DATA.Nx,DATA.Ny)))...
            .*real((ifft2(DATA.DENS_I(:,:,iz,it),DATA.Nx,DATA.Ny)));
        end
    end
    if REALP % plot in real space
        for it = FRAMES
            FIELD(:,:,it) = squeeze(compr(Gx(:,:,:,it)));
        end
    else % Plot spectrum
        process = @(x) abs(fftshift(x,2));
        shift_x = @(x) fftshift(x,2);
        shift_y = @(x) fftshift(x,2);
        tmp = zeros(DATA.Nx,DATA.Ny,Nz);
        for it = FRAMES
            for iz = 1:numel(DATA.z)
            tmp(:,:,iz) = process(squeeze(fft2(Gx(:,:,iz,it),DATA.Nx,DATA.Ny)));
            end
        FIELD(:,:,it) = squeeze(compr(tmp(1:DATA.Nkx,1:DATA.Nky,:)));
        end  
    end    
end
TOPLOT.FIELD     = FIELD;
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

