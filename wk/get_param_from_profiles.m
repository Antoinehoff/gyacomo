function [params, profiles] = get_param_from_profiles(folder,rho,Lref,mref,Bref,FROMPROFILE)
%get_param_from_profiles compute the input param for GYACOMO from profiles
%data
m_e           = 9.11e-31;
params.SIGMA  = sqrt(m_e/mref);
if ~FROMPROFILE
    nmvm = 5;
    stencil = 9;
    for field = {'ne','te','ti','wExB'}
        f_.name = field{1};
        fxy  = load([folder,f_.name,'.txt']);
        f_.x = movmean(fxy(:,1),nmvm); 
        f_.y = movmean(fxy(:,2),nmvm);
        f_.K =-calculate_derivative(f_.x,f_.y,stencil)./f_.y;
        profiles.(field{1}) = f_;
    end
    profiles.nmvm    = nmvm;
    profiles.stencil = stencil;
else
    profiles = read_DIIID_profile([folder,'/profiles.txt']);
end

    n_e  = interpolate_at_x0(profiles.ne.x,profiles.ne.y,rho);
    T_e  = interpolate_at_x0(profiles.te.x,profiles.te.y,rho);
    T_i  = interpolate_at_x0(profiles.ti.x,profiles.ti.y,rho);
    wExB = 0;%interpolate_at_x0(profiles.wExB.x,profiles.wExB.y,rho);
    K_w  = 0;%interpolate_at_x0(profiles.wExB.x,profiles.wExB.K,rho);
    params.K_Ne = interpolate_at_x0(profiles.ne.x,profiles.ne.K,rho);
    params.K_Ni = params.K_Ne;
    params.K_Te = interpolate_at_x0(profiles.te.x,profiles.te.K,rho);
    params.K_Ti = interpolate_at_x0(profiles.ti.x,profiles.ti.K,rho);
    params.TAU  = T_i/T_e;
    params.nuGENE = 2.3031E-5*Lref*(n_e)/(T_e)^2*(24.-log(sqrt(n_e*1.0E13)/T_e*0.001));
    params.NU =  3/8*sqrt(pi/2)*params.TAU.^(3/2)*params.nuGENE;
    params.BETA   = 403.0E-5*n_e*T_e/(Bref*Bref);
    cref= sqrt(T_e*1000*1.6e-19/mref);
    params.EXBRATE= K_w*wExB*1000*Lref/cref;
end