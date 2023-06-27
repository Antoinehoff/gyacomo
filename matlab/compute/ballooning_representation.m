function [field_chi,chi] = ballooning_representation(field,shear,Npol,kx,iky,ky,z)
    dims = size(field);
    Nkx  = dims(2);
    % Apply ballooning transform
    if(Nkx > 1)
        nexc = round(ky(2)*shear*2*pi/kx(2));
    else
        nexc = 1;
    end
    is   = max(1,iky-1);
    Npi  = (Nkx-1)-2*nexc*(is-1);
    if(Npi <= 1)
        ordered_ikx = 1;
    elseif(mod(Nkx,2) == 0)
        tmp_ = (Nkx-is+1):-is:(Nkx/2+2);
        ordered_ikx = [tmp_(end:-1:1), 1:is:(Nkx/2)];
    else
        Np_ = (Nkx+1)/(2*is);
        ordered_ikx = [(Np_+1):Nkx 1:Np_];
    end
    try
        idx=kx./kx(2);
    catch
        idx=0;
    end
    idx= idx(ordered_ikx);
    Nkx = numel(idx);

    field_chi = zeros(  Nkx*dims(3)  ,1);
    chi   = zeros(  Nkx*dims(3)  ,1);

    for i_ =1:Nkx
        start_ =  (i_-1)*dims(3) +1;
        end_ =  i_*dims(3);
        ikx  = ordered_ikx(i_);
        field_chi(start_:end_)  = field(iky,ikx,:);
    end

    % Define ballooning angle
    for i_ =1: Nkx
        for iz=1:dims(3)
            ii = dims(3)*(i_-1) + iz;
            chi(ii) =z(iz) + 2*pi*Npol*idx(i_)./is;
        end
    end
end
