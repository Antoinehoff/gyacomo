codir  = '/home/ahoffman/cosolver/';
matname= 'gk_landau_P16_J9_dk_5e-2_km_2.0_NFLR_8';
matdir = [codir,'LDGK_P16_J9_NFLR_8/matrices/'];
kperp  = load([matdir,'kperp.in']); 
Nkp    = numel(kperp);
outdir = '/home/ahoffman/gyacomo/iCa/';
outfile= [outdir,matname,'.h5'];

if(exist(outfile)>0)
    system(['rm ',outdir,matname,'.h5']);
end

Nmax = numel(kperp);
for n_=1:Nmax
    n_string = sprintf('%5.5d',n_); disp(n_string);
    filename  = ['ei.',n_string,'.h5'];
    if(exist([matdir,filename])>0)
        % Load matrices and write them in one h5 file

        olddname = '/Ceipj/CeipjF'; infile = [matdir,'ei.',n_string,'.h5'];
        newdname = ['/',sprintf('%5.5d',n_-1),olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname = '/Ceipj/CeipjT'; infile = [matdir,'ei.',n_string,'.h5'];
        newdname = ['/',sprintf('%5.5d',n_-1),olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname = '/Ciepj/CiepjT'; infile = [matdir,'ie.',n_string,'.h5'];
        newdname = ['/',sprintf('%5.5d',n_-1),olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname = '/Ciepj/CiepjF'; infile = [matdir,'ie.',n_string,'.h5'];
        newdname = ['/',sprintf('%5.5d',n_-1),olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname =  '/Caapj/Ceepj'; infile = [matdir,'self.',n_string,'.h5'];
        newdname = ['/',sprintf('%5.5d',n_-1),olddname];
        mat_ee   = load_write_mat(infile,olddname,outfile,newdname);

        olddname =  '/Caapj/Ciipj'; infile = [matdir,'self.',n_string,'.h5'];
        newdname = ['/',sprintf('%5.5d',n_-1),olddname];
        mat_ii   = load_write_mat(infile,olddname,outfile,newdname);

        % to verify symmetry
        verif = max(abs(imag(eig(mat_ee)))) + max(abs(imag(eig(mat_ii))));
        if(verif>0) disp(['Warning, non sym matrix at ', n_string]); end;
        % to verify negative EV
        verif = max(real(eig(mat_ee))) + max(real(eig(mat_ii)));
        if(verif>0) disp(['Warning, positive EV at ', n_string]); end;
    else
        break
    end
end
olddname = '/Ceipj/CeipjF'; infile = [matdir,'ei.',n_string,'.h5'];
Pmaxe    = h5readatt(infile,olddname,'Pmaxe');
Jmaxe    = h5readatt(infile,olddname,'Jmaxe');
Pmaxi    = h5readatt(infile,olddname,'Pmaxi');
Jmaxi    = h5readatt(infile,olddname,'Jmaxi');
dims_e   = [0 0];
dims_e(1)= Pmaxe; dims_e(2) = Jmaxe;
dims_i   = [Pmaxi; Jmaxi];
h5create(outfile,'/dims_e',numel(dims_e));
h5write (outfile,'/dims_e',dims_e);
h5create(outfile,'/dims_i',numel(dims_i));
h5write (outfile,'/dims_i',dims_i);
% 
h5create(outfile,'/coordkperp',numel(kperp(1:Nmax)));
h5write (outfile,'/coordkperp',kperp(1:Nmax)');
fid = H5F.open(outfile);

H5F.close(fid);
