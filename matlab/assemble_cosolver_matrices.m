codir  = '/home/ahoffman/cosolver/';
matname= 'gk.hacked_sugama_P_4_J_2_N_75_kpm_5.0';
matdir = [codir,matname,'/'];
kperp  = load([matdir,'kperp.log']); 
Nkp    = numel(kperp);
outdir = '/home/ahoffman/HeLaZ/iCa/';
outfile= [outdir,matname,'.h5'];

if(exist(outfile)>0)
    system(['rm ',outdir,matname,'.h5']);
end

Nmax = numel(kperp);
for n_=1:Nmax
    n_string = sprintf('%5.5d',n_-1); disp(n_string);
    scandir  = ['scanfiles_',n_string,'/'];
    if(exist([matdir,scandir,  'ei.h5'])>0)
        % Load matrices and write them in one h5 file

        olddname = '/Ceipj/CeipjF'; infile = [matdir,scandir,  'ei.h5'];
        newdname = ['/',n_string,olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname = '/Ceipj/CeipjT'; infile = [matdir,scandir,  'ei.h5'];
        newdname = ['/',n_string,olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname = '/Ciepj/CiepjT'; infile = [matdir,scandir,  'ie.h5'];
        newdname = ['/',n_string,olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname = '/Ciepj/CiepjF'; infile = [matdir,scandir,  'ie.h5'];
        newdname = ['/',n_string,olddname];
        load_write_mat(infile,olddname,outfile,newdname);

        olddname =  '/Caapj/Ceepj'; infile = [matdir,scandir,'self.h5'];
        newdname = ['/',n_string,olddname];
        mat_ee   = load_write_mat(infile,olddname,outfile,newdname);

        olddname =  '/Caapj/Ciipj'; infile = [matdir,scandir,'self.h5'];
        newdname = ['/',n_string,olddname];
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
olddname = '/Ceipj/CeipjF'; infile = [matdir,scandir,  'ei.h5'];
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