Nkx   = 6;
Nky   = 3;

my    = 0:(Nky-1);
mx    = zeros(1,Nkx);

PERIODIC = 1;
Npol  = 1;
Nexc  = 0;

shear = 0.8;
Ly    = 120;
Lx    = 120;
dky   = 2*pi/Ly;


for ix = 1:Nkx
    if(mod(Nkx,2) == 0)%even
        mx_max  = (Nkx/2);
        if(ix<=Nkx/2+1)
            mx(ix) = (ix-1);
        else
            mx(ix) = ix-Nkx-1;
        end
    else %odd
        mx_max  = (Nkx-1)/2;
        if(ix<=(Nkx-1)/2+1)
            mx(ix) = (ix-1);
        else
            mx(ix) = ix-Nkx-1;
        end        
    end
end
disp(mx)


if Nexc == 0 %% Adapt Nexc
    dkx = 2*pi/Lx;
    Nexc = ceil(0.9*2*pi*shear*dky/dkx);
else
    dkx   = 2*pi*shear*dky/Nexc;
end

kx = mx*dkx;
ky = my*dky;

kx_max = mx_max*dkx;
ikx_zBC_R = zeros(Nky,Nkx);
for iy = 1:Nky
    shift = 2*pi*shear*ky(iy)*Npol;
    for ix = 1:Nkx
        kx_shift = kx(ix) + shift;
        if ((kx_shift > kx_max) && (~PERIODIC))
            ikx_zBC_R(iy,ix) = nan;
        else
            ikx_zBC_R(iy,ix) = ix+(iy-1)*Nexc;
         if(ikx_zBC_R(iy,ix) > Nkx)
%             ikx_zBC_R(iy,ix) = ikx_zBC_R(iy,ix) - Nkx;
            ikx_zBC_R(iy,ix) = mod(ikx_zBC_R(iy,ix)-1,Nkx)+1;
         end
        end
    end 
end
disp(ikx_zBC_R)

kx_min = (-mx_max+(1-mod(Nkx,2)))*dkx;
ikx_zBC_L = zeros(Nky,Nkx);
for iy = 1:Nky
    shift = 2*pi*shear*ky(iy)*Npol;
    for ix = 1:Nkx
        kx_shift = kx(ix) - shift;
        if ((kx_shift < kx_min) && (~PERIODIC))
            ikx_zBC_L(iy,ix) = nan;
        else
            ikx_zBC_L(iy,ix) = ix-(iy-1)*Nexc;
         if(ikx_zBC_L(iy,ix) < 1)
%             ikx_zBC_L(iy,ix) = ikx_zBC_L(iy,ix) + Nkx;
            ikx_zBC_L(iy,ix) = mod(ikx_zBC_L(iy,ix)-1,Nkx)+1;
         end
        end
    end 
end
disp(ikx_zBC_L)

