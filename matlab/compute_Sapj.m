function [ Sapj ] = compute_Sapj(p, j, Kr, Kz, Napj, specie, phi, MODEL, GRID)
%COMPUT_SAPJ compute the non linear term for moment pj specie a
%   perform a convolution by product of inverse fourier transform and then
%   put the results back into k space with an FFT
    Jmax = GRID.jmaxe * (specie=='e')...
         + GRID.jmaxi * (specie=='i');
    %padding
    Pad = 2.0;
    Sapj = zeros(numel(Kr), numel(Kz));
    F    = zeros(numel(Kr), numel(Kz));
    G    = zeros(numel(Kr), numel(Kz));
    for n = 0:Jmax % Sum over Laguerre

        for ikr = 1:numel(Kr)
            for ikz = 1:numel(Kz)
                kr = Kr(ikr); kz = Kz(ikz); 
                BA = sqrt(kr^2+kz^2)*...
                      (MODEL.sigma_e*sqrt(2) * (specie=='e')...
                      + sqrt(2*MODEL.tau_i)  * (specie=='i'));
                  
                %First conv term
                F(ikr,ikz) = (kz-kr).*phi(ikr,ikz).*kernel(n,BA);
                %Second conv term
                G(ikr,ikz) = 0.0;
                for s = 0:min(n+j,Jmax)
                    G(ikr,ikz) = ...
                        G(ikr,ikz) + dnjs(n,j,s) .* squeeze(Napj(p+1,s+1,ikr,ikz));
                end
                G(ikr,ikz) = (kz-kr) .* G(ikr,ikz);
            end
        end
        %Conv theorem
        Sapj = Sapj + conv_thm_2D(F,G,Pad);
    end
end


