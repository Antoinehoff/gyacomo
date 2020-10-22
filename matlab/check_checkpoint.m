function [] = check_checkpoint(cp_ID, OUTPUTS)

cp_fname = OUTPUTS.rstfile0;
cp_fname = [cp_fname(2:end-1),'_%.2d.h5'];
cp_fname = sprintf(cp_fname,cp_ID)

tmp      = h5read(cp_fname,'/Basic/moments_i');
Nipj_cp  = tmp.real + 1i * tmp.imaginary;
t_cp     = h5readatt(cp_fname,'/Basic','time');

[~,~,Nr,Nz] = size(Nipj_cp);  

% ni00_cp = fftshift(ifft2(half_2_full_cc_2D(squeeze(Nipj_cp(1,1,:,:))),'symmetric'));
ni00_cp = real(fftshift(ifft2(squeeze(Nipj_cp(1,1,:,:)),Nz,Nz)));
fig  = figure;
    pclr=pcolor(transpose(ni00_cp(:,:)));set(pclr, 'edgecolor','none'); colorbar;
    title(['$n_i^{00}$',', $t \approx$', sprintf('%.3d',ceil(t_cp))]);
    xlabel('$i_r$'); ylabel('$i_z$');
end