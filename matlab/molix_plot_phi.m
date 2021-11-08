function [phib, b_angle] = molix_plot_phi(resfile,varargin)
% plot ballooning representation at the specified iput by time2plot for
% each ky modes


set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16);
set(0, 'DefaultLineLineWidth', 1.5);


% read data and attributes
coordkx = h5read(resfile,'/data/var2d/phi/coordkx');
dimskx = size(coordkx);
Nkx = (dimskx(1) - 1)/2;
% Nkx = (length(coordkx)-1)/2
coordz = h5read(resfile,'/data/var2d/phi/coordz');
coordky = h5read(resfile,'/data/var2d/phi/coordky');
Nky = length(coordky);
Nz = length(coordz);
coordtime =h5read(resfile,'/data/var2d/time');
Nt = length(coordtime);
pmax = double(h5readatt(resfile,'/data/input','pmaxi'));
jmax = double(h5readatt(resfile,'/data/input','jmaxi'));
epsilon = double(h5readatt(resfile,'/data/input','inv. aspect ratio'));
safety_fac = double(h5readatt(resfile,'/data/input','safety factor'));
nu  = double(h5readatt(resfile,'/data/input','nu'));
shear = double(h5readatt(resfile,'/data/input','magnetic shear'));

RTi= double(h5readatt(resfile,'/data/input','RTi'));
RTe= double(h5readatt(resfile,'/data/input','RTe'));
Rn= double(h5readatt(resfile,'/data/input','Rni'));

% read electrostatic pot.
if(nargin == 1 ) 
      % plot end of the simulation
     if( Nt >1)
         [~,idxte] = min(abs(coordtime -coordtime(end-1)));
     else
        idxte =0;
     end
elseif(nargin ==2)
    time2plot = varargin{1};
    [~,idxte] = min(abs(coordtime -time2plot));
end

iframe = idxte-1;

dataset = ['/data/var2d/phi/',num2str(iframe,'%06d')];
phi = h5read(resfile,dataset);

% read eigenvalu
% dataset = ['/data/var2d/eig/growth_rate'];
% growth_rate = h5read(resfile,dataset);
% dataset = ['/data/var2d/eig/freq'];
% frequency =  h5read(resfile,dataset);

dataset = ['/data/var2d/time'];
time2d= h5read(resfile,dataset);

% Apply baollooning tranform
for iky=1:Nky
    dims = size(phi.real);
    phib_real = zeros(  dims(1)*Nz  ,1);
    phib_imag= phib_real;
    b_angle = phib_real;
    
    midpoint = floor((dims(1)*Nz )/2)+1;
    for ip =1: dims(1)
        start_ =  (ip -1)*Nz +1;
        end_ =  ip*Nz;
        phib_real(start_:end_)  = phi.real(ip,iky,:);
        phib_imag(start_:end_)  = phi.imaginary(ip,iky, :);
    end
    
    % Define ballooning angle
    idx = -Nkx:1:Nkx;
    for ip =1: dims(1)
        for iz=1:Nz
            ii = Nz*(ip -1) + iz;
            b_angle(ii) =coordz(iz) + 2*pi*idx(ip);
        end
    end
    
    % normalize real and imaginary parts at chi =0
    [~,idxLFS] = min(abs(b_angle -0));
    phib = phib_real(:) + 1i * phib_imag(:);
    % Normalize to the outermid plane
    phib_norm = phib(:);% / phib( idxLFS)    ;
    phib_real_norm(:)  = real( phib_norm);%phib_real(:)/phib_real(idxLFS);
    phib_imag_norm(:)  = imag( phib_norm);%phib_imag(:)/ phib_imag(idxLFS);
%     
%     
%     figure; hold all;
%     plot(b_angle/pi, phib_real_norm,'-b');
%     plot(b_angle/pi, phib_imag_norm ,'-r');
%     plot(b_angle/pi, sqrt(phib_real_norm .^2 + phib_imag_norm.^2),'--k');
%     legend('real','imag','norm')
%     xlabel('$\chi / \pi$')
%     ylabel('$\phi_B (\chi)$');
% %     title(['$(P,J) =(',num2str(pmax),', ', num2str(jmax),'$), $\nu =',num2str(nu),...
% %         '$, $\epsilon = ',num2str(epsilon),'$, $k_y = ', num2str(coordky( iky)),'$, $q =',num2str(safety_fac),'$, $s = ', num2str(shear),'$, $R_N = ', ...
% %         num2str(Rn),'$, $R_{T_i} = ', num2str(RTi),'$, $N_z =',num2str(Nz),'$']);
%     title(['molix, t=',num2str(coordtime(idxte))]);
%     %set(gca,'Yscale','log')
    %
    
%     %Check symmetry of the mode at the outter mid plane
%     figure; hold all;
%     right = phib_real(midpoint+1:end);
%     left = fliplr(phib_real(1:midpoint-1)');
%     up_down_symm  = right(1:end) - left(1:end-1)';
%     %plot(b_angle(midpoint+1:end)/pi,phib_real(midpoint+1:end),'-xb');
%     plot(b_angle(midpoint+1:end)/pi,up_down_symm ,'-xb');
    %plot(abs(b_angle(1:midpoint-1))/pi,phib_real(1:midpoint-1),'-xb');
    %
    %
    % figure; hold all
    % plot(b_angle/pi, phib_imag.^2 + phib_real.^2 ,'xk');
    % %set(gca,'Yscale','log')
    % xlabel('$\chi / \pi$')
    % ylabel('$\phi_B (\chi)$');
    % title(['$(P,J) =(',num2str(pmax),', ', num2str(jmax),'$), $\nu =',num2str(nu),'$, $\epsilon = ',num2str(epsilon),'$, $q =',num2str(safety_fac),'$, $s = ', num2str(shear),'$, $k_y =',num2str(ky),'$']);
end

end