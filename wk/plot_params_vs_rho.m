function [ ] = plot_params_vs_rho(geom,prof_folder,rho,rho_min,rho_max,Lref,mref,Bref,FROMPROFILE)

if FROMPROFILE
    % profiles = read_DIIID_profile([prof_folder,'/profiles.txt']);
    [params,profiles] = get_param_from_profiles(prof_folder,rho,Lref,mref,Bref,FROMPROFILE);
    m_e    = 9.11e-31;
    sigma  = sqrt(m_e/mref);
    TAU_a  = profiles.ti.y./profiles.te.y;
    nuGENE = 2.3031E-5*Lref*(profiles.ne.y)./(profiles.te.y).^2 ...
    .*(24.-log(sqrt(profiles.ne.y*1.0E13)./profiles.te.y*0.001));
    NU_a   = 3/8*sqrt(pi/2)*TAU_a.^(3/2).*nuGENE;
    BETA_a = 403.0E-5*profiles.ne.y.*profiles.te.y/(Bref*Bref);
    figure
    subplot(231)
        plot(profiles.ne.x,profiles.ne.y); hold on
        xlabel('$\rho$'); ylabel('$n_e(\times 10^{19}m^{-1})$')
        xlim([rho_min rho_max]); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        grid on
    subplot(234)
        plot(profiles.te.x,profiles.te.y,'DisplayName','electrons'); hold on
        plot(profiles.ti.x,profiles.ti.y,'DisplayName','ions'); hold on
        xlabel('$\rho$'); ylabel('$T$(keV)')
        xlim([rho_min rho_max]); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        grid on
    subplot(235)
    yyaxis("left")
        semilogy(profiles.ti.x,NU_a); hold on
        ylabel('$\nu$')
    yyaxis("right")
        plot(profiles.ti.x,BETA_a*100); hold on
        ylabel('$\beta$[\%]')
        xlim([rho_min rho_max]);xlabel('$\rho$'); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        grid on
    subplot(232)
        plot(profiles.te.x,profiles.te.K,'DisplayName','$\kappa_{Te}$'); hold on
        plot(profiles.ti.x,profiles.ti.K,'DisplayName','$\kappa_{Ti}$'); hold on
        plot(profiles.ne.x,profiles.ne.K,'DisplayName','$\kappa_{n}$'); hold on
        xlim([rho_min rho_max]); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        xlabel('$\rho$'); ylabel('$\kappa_\chi=R_0/L_\chi$')
        legend show
        grid on
else
    rho_a = linspace(rho_min,rho_max,100);
    
    NU_a   = zeros(size(rho_a));
    BETA_a = zeros(size(rho_a));
    K_Ne_a = zeros(size(rho_a));
    K_Ti_a = zeros(size(rho_a));
    K_Te_a = zeros(size(rho_a));
    
    for i = 1:numel(rho_a)
        rho_ = rho_a(i);
        [params,profiles] = get_param_from_profiles(prof_folder,rho_,Lref,mref,Bref,FROMPROFILE);
        NU_a(i)   = 3/8*sqrt(pi/2)*params.TAU.^(3/2)*params.nuGENE;
        BETA_a(i) = params.BETA;
        K_Ne_a(i) = params.K_Ne;
        K_Ti_a(i) = params.K_Ti;
        K_Te_a(i) = params.K_Te;
    end
    figure
    subplot(231)
        plot(profiles.ne.x,profiles.ne.y); hold on
        xlabel('$\rho$'); ylabel('$n_e(\times 10^{19}m^{-1})$')
        xlim([rho_min rho_max]); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        grid on
    subplot(234)
        plot(profiles.te.x,profiles.te.y,'DisplayName','electrons'); hold on
        plot(profiles.ti.x,profiles.ti.y,'DisplayName','ions'); hold on
        xlabel('$\rho$'); ylabel('$T$(keV)')
        xlim([rho_min rho_max]); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        grid on
    subplot(235)
    yyaxis("left")
        semilogy(rho_a,NU_a); hold on
        ylabel('$\nu$')
    yyaxis("right")
        plot(rho_a,BETA_a*100); hold on
        ylabel('$\beta$[\%]')
        xlim([rho_min rho_max]);xlabel('$\rho$'); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        grid on
    subplot(232)
        plot(rho_a,K_Te_a,'DisplayName','$\kappa_{Te}$'); hold on
        plot(rho_a,K_Ti_a,'DisplayName','$\kappa_{Ti}$'); hold on
        plot(rho_a,K_Ne_a,'DisplayName','$\kappa_{n}$'); hold on
        xlim([rho_min rho_max]); xline(rho,'--k','Linewidth',2,'DisplayName','$\rho_{sim}$')
        xlabel('$\rho$'); ylabel('$\kappa_\chi=R_0/L_\chi$')
        legend show
        grid on
end
    subplot(133)
    [R,Z] = plot_miller(geom,rho,32,0);
    plot(R,Z,'-b');
    xlabel('R [m]');
    ylabel('Z [m]');
    axis tight
    axis equal
end