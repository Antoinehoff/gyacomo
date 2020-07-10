%meshgrid
[KR, KZ] = meshgrid(kr,kz);
%BE = sqrt(KR.^2+KZ.^2)*params.mu*sqrt(2);
BI = sqrt(KR.^2+KZ.^2)*sqrt(2*params.tau);

%non linear term, time evolution
Sipj     = zeros(GRID.nkr, GRID.nkz, numel(time));

%padding
Mr = 2*GRID.nkr; Mz = 2*GRID.nkz;
F_pad    = zeros(Mr, Mz);
G_pad    = zeros(Mr, Mz);

p = 0; j = 0; %pj moment

for it = 1:numel(time) % time loop
    Sipj_pad = zeros(Mr, Mz);
    for n = 0, GRID.jmaxi % Sum over Laguerre
        %First conv term
        F = (KZ-KR).*phiHeLaZ(it,:,:)*kernel(n,BE);
        
        %Second conv term
        G = 0.*(KZ-KR);
        for s = 0:min(n+j,GRID.jmaxi)
            G = G + dnjs(n,j,s) .* Napj(p,s,:,:,it);
        end
        G = (KZ-KR) .* G;
        
        %Padding
        F_pad(1:GRID.nkr,1:GRID.nkz) = F_;
        G_pad(1:GRID.nkr,1:GRID.nkz) = G_;
        
        %Inverse fourier transform to real space
        f = ifft2(F_pad);
        g = ifft2(G_pad);
        
        %Conv theorem
        Sipj_pad = Sipj_pad + fft2(f.*g);
    end
    Sipj(:,:,it) = Sipj_pad(1:GRID.nkr,1:GRID.nkz);
end