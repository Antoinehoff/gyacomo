% Hflux_x = 0;
% Hflux_x = 0 * data.Ts5D;
filename = '/home/ahoffman/HeLaZ/results/shearless_cyclone/64x32x16x5x3_CBC_100/outputs_00.h5';
kernel_i = h5read(filename,'/data/metric/kernel_i');
Jacobian = h5read(filename,'/data/metric/Jacobian');
STEPS = 1:numel(data.Ts5D);
Hflux_x = 1:numel(STEPS);
its_ = 1;
for it = STEPS
    t = data.Ts5D(it);
    [~,it5d] = min(abs(data.Ts5D-t));
    [~,it3d] = min(abs(data.Ts3D-t));

    Nj = data.Jmaxi;
    Nz = data.Nz; z  = data.z; dz = data.z(2)-data.z(1);
    kx = data.kx; ky = data.ky; Nkx = data.Nkx; Nky = data.Nky;
    Jz = squeeze(Jacobian(:,1));
    invjac = 1/2/pi/data.Q0;
    
    % Factors for sum kern mom
    c2n   = @(n_) 0.5*sqrt(2);
    c0n   = @(n_) 2*n_+ 1.5;
%     c0n   = @(n_) 2*n_+ 2/3;
    c0np1 = @(n_) -(n_+1);
    c0nm1 = @(n_) -n_;
    
    % Factors for correction operator
%     dn   = @(n_) -2*(n_+ 1.5);
    dn   = @(n_) -2*n_+ 1.5;
    dnp1 = @(n_) (n_+1);
    dnm1 = @(n_) n_;
    
    % BUILD TERM TO SUM
    sumkernmom = zeros(Nky,Nkx,Nz);
    correct_op = zeros(Nky,Nkx,Nz);
    for in = 1:Nj
        n = in-1;
        Kn    = squeeze(kernel_i(in,:,:,:,1));
        N2n   = squeeze(data.Nipj(3,in  ,:,:,:,it5d));
        N0n   = squeeze(data.Nipj(1,in  ,:,:,:,it5d));
        sumkernmom = sumkernmom + ...
            Kn.* (c2n(n) .* N2n + c0n(n) .* N0n);
        
        correct_op = correct_op + ...
            Kn.* dn(n) .* Kn;
        
        if(in > 1)
            N0nm1 = squeeze(data.Nipj(1,in-1,:,:,:,it5d));
            sumkernmom = sumkernmom + Kn.* c0nm1(n).*N0nm1;
            
            Knm1  = squeeze(kernel_i(in-1,:,:,:,1));
            correct_op = correct_op + Kn.* dnm1(n) .* Knm1;       
        end
        if(in<Nj)
            N0np1 = squeeze(data.Nipj(1,in+1,:,:,:,it5d));
            sumkernmom = sumkernmom + Kn.* c0np1(n).*N0np1;
            
            Knp1  = squeeze(kernel_i(in+1,:,:,:,1));
            correct_op = correct_op + Kn.* dnp1(n) .* Knp1;       
        end
    end
    [~,KYY,ZZZ] =  meshgrid(kx,ky,z);
    % -- adding correction term
    correction_term = data.PHI(:,:,:,it3d) .* correct_op/ data.scale;
    
    % -- summing up
    u =  -1i*KYY.*data.PHI(:,:,:,it3d);
    n =   sumkernmom + correction_term;
    clear sumkernmom correction_term correct_op KYY;
    
    sum_kxky  = conj(u).*n;
    half_plane= sum_kxky;
    Ny = 2*Nky-1;
    sum_kxky  = zeros(Ny,Nkx,Nz);
    sum_kxky(1:Nky,:,:)=half_plane(:,:,:);
    sum_kxky((Nky+1):(Ny),1,:)=conj(half_plane(Nky:-1:2,1,:));
    sum_kxky((Nky+1):(Ny),2:Nkx,:)=conj(half_plane(Nky:-1:2,Nkx:-1:2,:));  
    
    sum_kxky  = squeeze(sum(sum(sum_kxky,1),2));
    % Z average
    J_= 0;
    q_ = 0;   
    for iz = 1:Nz
%         add = u(1,1,iz)*n(1,1,iz)...
%           +2*sum(real(u(2:Nky,1,iz).*conj(n(2:Nky,1,iz))),1)...
%           +sum(2*real(u(1,2:Nkx/2+1,iz).*conj(n(1,2:Nkx/2+1,iz))) ...
%             +sum(2*real(u(2:Nky,2:Nkx/2+1,iz).*conj(n(2:Nky,2:Nkx/2+1,iz)))+...
%                  2*real(u(2:Nky,Nkx/2+1:Nkx,iz).*conj(n(2:Nky,Nkx/2+1:Nkx,iz))),1),2);
        add = sum_kxky(iz);
        if(mod(iz,2)==1)
            q_ = q_ + 4 * Jz(iz)*add;
            J_ = J_ + 4 * Jz(iz);
        else
            q_ = q_ + 2 * Jz(iz)*add;
            J_ = J_ + 2 * Jz(iz);
        end
    end
    q_ = q_/J_;
    Hflux_x(its_) = real(q_);
    its_ = its_ + 1;
end
%%
figure
plot(data.Ts0D,data.HFLUX_X.*data.scale); hold on;
plot(data.Ts5D(STEPS),Hflux_x.*data.scale)