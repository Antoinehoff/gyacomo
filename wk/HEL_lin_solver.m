% We formulate the linear system as dN/dt + A*N = 0

q  = 1;
Jxyz = 1; hatB = 1; dzlnB = 0;
ddz  = 0;
tau = 0.001;
IVAN = 1;
%ordering
O1_n    = 1;
O1_u    = 1-IVAN;
O1_Tpar = 1-IVAN;
O1_Tper = 1-IVAN;

RN = 0;
RT = 1*2*q/tau;
NU = 0.1*3/2/tau/(4-IVAN*2.25);
MU = 0.0;

kx_a = 0;linspace(-2.0,2.0,256);
ky_a = linspace(0.01,3.5,256);
g_a  = zeros(numel(kx_a),numel(ky_a)); w_a = g_a;
for i = 1:numel(kx_a)
for j = 1:numel(ky_a)
    kx    = kx_a(i);
    ky    = ky_a(j);
    lperp = (kx*kx + ky*ky)/2;
    Cperp =-1i*ky;
    if IVAN
        K0    = 1 - lperp*tau;% + 1/2*lperp^2*tau^2;
        K1    = lperp*tau;% - lperp^2*tau^2;
    else
        K0    = exp(-lperp*tau);
        K1    = lperp*tau*exp(-lperp*tau);
    end
    
    % Magnetic terms
    M_n_n    = 2*tau*Cperp/q*O1_n;
    M_n_u    = 0;
    M_n_Tpar = sqrt(2)*tau*Cperp/q*O1_n;
    M_n_Tper =-tau*Cperp/q*O1_n;
    M_n_phi  = (2*K0-K1*O1_n)*Cperp;

    M_u_n    = 0;
    M_u_u    = 4*tau*Cperp/q*O1_u;
    M_u_Tpar = 0;
    M_u_Tper = 0;
    M_u_phi  = 0;

    M_Tpar_n    = sqrt(2)*tau*Cperp/q*O1_Tpar;
    M_Tpar_u    = 0;
    M_Tpar_Tpar = 6*tau*Cperp/q*O1_Tpar;
    M_Tpar_Tper = 0;
    M_Tpar_phi  = sqrt(2)*K0*Cperp;

    M_Tper_n    = -tau*Cperp/q*O1_Tper;
    M_Tper_u    = 0;
    M_Tper_Tpar = 0;
    M_Tper_Tper = 4*tau*Cperp/q*O1_Tper;
    M_Tper_phi  = -(K0-4*K1)*tau*Cperp*O1_Tper;

    M = [M_n_n,    M_n_u,    M_n_Tpar,    M_n_Tper,    M_n_phi;...
         M_u_n,    M_u_u,    M_u_Tpar,    M_u_Tper,    M_u_phi;...
         M_Tpar_n, M_Tpar_u, M_Tpar_Tpar, M_Tpar_Tper, M_Tpar_phi;...
         M_Tper_n, M_Tper_u, M_Tper_Tpar, M_Tper_Tper, M_Tper_phi];

    % Diamagnetic (gradients)
    D_n_n    = 0;
    D_n_u    = 0;
    D_n_Tpar = 0;
    D_n_Tper = 0;
    D_n_phi  = 1i*ky*(K0*RN - K1*RT*O1_n);

    D_u_n    = 0;
    D_u_u    = 0;
    D_u_Tpar = 0;
    D_u_Tper = 0;
    D_u_phi  = 0;

    D_Tpar_n    = 0;
    D_Tpar_u    = 0;
    D_Tpar_Tpar = 0;
    D_Tpar_Tper = 0;
    D_Tpar_phi  = 1i*ky*K0*RT/sqrt(2);

    D_Tper_n    = 0;
    D_Tper_u    = 0;
    D_Tper_Tpar = 0;
    D_Tper_Tper = 0;
    D_Tper_phi  = 1i*ky*(-K0*RT + K1*(RN+2*RT)*O1_Tper);

    D = [D_n_n,    D_n_u,    D_n_Tpar,    D_n_Tper,    D_n_phi;...
         D_u_n,    D_u_u,    D_u_Tpar,    D_u_Tper,    D_u_phi;...
         D_Tpar_n, D_Tpar_u, D_Tpar_Tpar, D_Tpar_Tper, D_Tpar_phi;...
         D_Tper_n, D_Tper_u, D_Tper_Tpar, D_Tper_Tper, D_Tper_phi];

    % Collision (GK Lenhard-Bernstein)
    CLB      = zeros(4,5);
    CLB(1,1) =-2*tau*lperp*O1_n;     
    CLB(1,5) =-2*q*K0*lperp;
    CLB(2,2) =-(1+2*tau*lperp*O1_u);
    CLB(3,3) =-(2+2*tau*lperp*O1_Tpar);
    CLB(4,4) =-(2+2*tau*lperp*O1_Tper); 
    CLB(4,5) =-q*K1*(2+2*tau*lperp)/tau*O1_Tper;
    % Collision (Dougherty terms)
    CDG      = zeros(4,5);
    CDG(1,1) = 4/3*K1^2+2*tau*K0^2*lperp*O1_n;
    CDG(1,3) =-2/3*sqrt(2)*K0*K1*O1_n;
    CDG(1,4) = 4/3*(K0-2*K1)*K1*O1_n + 2*tau*K0*(K1-K0)*lperp*O1_n;
    CDG(1,5) = 8/3/tau*(K0-K1)*K1^2*O1_n + 2*q*K0*(K0^2-K0*K1+K1^2)*lperp;
    CDG(2,2) = K0^2;
    CDG(3,1) =-2/3*sqrt(2)*K0*K1*O1_Tpar;
    CDG(3,3) = 2/3*K0^2;
    CDG(3,4) =-2/3*sqrt(2)*K0*(K0-2*K1*O1_Tpar);
    CDG(3,5) = 4/3/tau*sqrt(2)*K0*K1*(K1*O1_Tpar-K0);
    CDG(4,1) = 4/3*(K0-2*K1)*K1*O1_Tper + 2*tau*K0*(K1*O1_Tper-K0)*lperp;
    CDG(4,3) =-2/3*sqrt(2)*K0*(K0-2*K1*O1_Tper);
    CDG(4,4) = 4/3*(K0-2*K1*O1_Tper)^2 + 2*tau*(K0-K1)^2*lperp*O1_Tper;
    CDG(4,5) =-2/3*q/tau*(K0-K1)*(3*tau*K0^2*lperp-K0*K1*(4+3*tau*lperp)+K1^2*(8+3*tau*lperp));
    if 1
        C = NU*(CLB + CDG);
    else %to check
        C = zeros(4,5);
        C(1,1) = -2/3*tau*lperp*(4*tau*lperp);
        C(1,3) = -2/3*tau*lperp*(sqrt(2)-2*sqrt(2)*tau*lperp);
        C(1,4) = -2/3*tau*lperp*(1-tau*lperp);
        C(1,5) = -2/3*tau*lperp*(5*q*lperp - 11*q*tau*lperp^2);
        C(2,2) = -(1+2*tau*lperp);
        C(3,1) = -2/3*(sqrt(2)*tau*lperp);
        C(3,3) = -2/3*(2+5*tau*lperp);
        C(3,4) = -2/3*(sqrt(2)-4*sqrt(2)*tau*lperp);
        C(3,5) = -2/3*(2*sqrt(2)*q*lperp+8*sqrt(2)*q*tau*lperp^2);
        C(4,1) = -2/3*(tau*lperp - tau^2*lperp^2);
        C(4,3) = -2/3*(sqrt(2)*4*sqrt(2)*tau*lperp);
        C(4,4) = -2/3*(1+12*tau*lperp);
        C(4,5) = -2/3*(2*q*lperp+9*q*tau*lperp^2);
        C      = NU*C;

    end
    % Numerical diffusion
    Diff = MU*lperp^2*eye(4);

    % Build the matrix
    A = zeros(4);
    if 1
        A = M(1:4,1:4) + D(1:4,1:4) - C(1:4,1:4) + Diff;
        P = M(:,5) + D(:,5) - C(:,5);
    else
        A(1,1:4) = [      2*tau*Cperp/q,             0, sqrt(2)*tau*Cperp/q,  -tau*Cperp/q];
        A(2,1:4) = [                  0, 4*tau/q*Cperp,                   0,             0];
        A(3,1:4) = [sqrt(2)*tau*Cperp/q,             0,       6*tau*Cperp/q,             0];
        A(4,1:4) = [       -tau*Cperp/q,             0,                   0, 4*tau*Cperp/q];
        P = [ 1i*ky*(K0*RN-K1*RT) + (2*K0-K1)*Cperp;...
                                                  0; ...
                    K0*(1i*ky*RT + 2*Cperp)/sqrt(2); ...
              1i*ky*(-K0*RT+K1*(RN+2*RT))-(K0-4*K1)*Cperp];
    end
    % Solve Poisson (op_phi * phi = charge_i + charge_e
    op_phi   = q^2/tau*(1-(K0^2 + K1^2)) + 1;

    A(:,1) = A(:,1) + q*K0/op_phi.*P;
    A(:,4) = A(:,4) + q*K1/op_phi.*P;

    % the system is written as -iwf + Af = 0 <-> iwf = lambda f where
    % lambda is an eigenvalue of A
    lambda = eig(A);
    gamma  =-real(lambda);
    omega  = imag(lambda);
    [g_a(i,j),idx] = max(gamma);
    % w_a(i,j) = omega(idx);
    w_a(i,j) = max(omega);
end
end
%
figure
if (numel(kx_a) == 1)
    plot(ky_a,g_a); hold on
    % plot(ky_a,w_a,'--'); hold on
    xlabel('$k_y\rho_s$',Interpreter='latex',FontSize=18);
    ylabel('$\gamma R_0/c_s$',Interpreter='latex',FontSize=18);
    fontsize(18,"points")
else
    contourf(kx_a,ky_a,g_a',50)
    ylabel('$k_y\rho_s$',Interpreter='latex',FontSize=18);
    xlabel('$k_x\rho_s$',Interpreter='latex',FontSize=18);
    colorbar
    clim(0.4*[-1 1])
    colormap(bluewhitered)
end