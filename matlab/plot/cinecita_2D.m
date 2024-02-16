function [MOVIE] = cinecita_2D(DATADIR,fieldname,varargin)
DATADIR = [DATADIR,'/'];
J0 = 0; J1 = 20; nskip = 1; FOURIER = 0; SAT = 0.75;
PLAY = 1;
if numel(varargin) > 0
    if isfield(varargin{1},'J0') 
        J0 = varargin{1}.J0;
    end
    if isfield(varargin{1},'J1') 
        J1 = varargin{1}.J1;
    end
    if isfield(varargin{1},'nskip') 
        nskip = varargin{1}.nskip;
    end
    if isfield(varargin{1},'FOURIER') 
        FOURIER = varargin{1}.FOURIER;
    end
    if isfield(varargin{1},'SAT') 
        SAT = varargin{1}.SAT;
    end
    if isfield(varargin{1},'PLAY') 
        PLAY = varargin{1}.PLAY;
    end
end
%
data = {};
data = compile_results_low_mem(data,DATADIR,J0,J1);
data = load_params(data,data.outfilenames{1});
Nx = data.grids.Nx; Ny = data.grids.Ny; Na = data.inputs.Na;
[KX,KY] = meshgrid(data.grids.kx,data.grids.ky);
format = '.gif';
%% Select the field
OPE_         = 1; % Operation on data
Nspec        = 0; % Load a specie
switch fieldname
case '\phi' %ES pot
    toload  = 'phi';
    ltxname = '\phi';
    vname   = 'phi';
case '\phi_{zf}' %ES pot
    toload  = 'phi';
    OPE_    = KY==0; % Operation on data
    ltxname = '\phi_{zf}';    
    vname   = 'phizf';
case '\psi' %EM pot
    toload  = 'psi';
    ltxname = '\psi';
    vname   = 'psi';
case 'n_i' % ion density
    toload  = 'dens';
    Nspec   = 1;
    ltxname = 'n_i';
    vname   = 'dens_i';
case 'n_e' % ion density
    toload  = 'dens';
    Nspec   = 2;
    ltxname = 'n_e';
    vname   = 'dens_e';
case 'v_{Ex}' %ES pot
    toload = 'phi';
    OPE_   = -1i*KY; % Operation on data
    ltxname= 'v_{Ex}';
    vname   = 'vEx';
case 'v_{Ey}' %ES pot
    toload = 'phi';
    OPE_   = 1i*KX; % Operation on data
    ltxname= 'v_{Ey}';
    vname   = 'vEy';
case '\delta B_x'
    toload = 'psi';
    OPE_   = -1i*KY; % Operation on data
    ltxname = '\delta B_x';
    vname   = 'dBx';
case '\delta B_y'
    toload = 'psi';
    OPE_   = 1i*KX; % Operation on data
    ltxname = '\delta B_y';    
    vname   = 'dBy';
otherwise
    toload = fieldname;
    OPE_   = 1;
end

if Nspec
    [FIELD,TIME] = compile_results_2Da(DATADIR,J0,J1,[toload,'_obmp']);
    FIELD = squeeze(FIELD(Nspec,:,:,:));
else
    [FIELD,TIME] = compile_results_2D(DATADIR,J0,J1,[toload,'_obmp']);
end

TIME = TIME(1:nskip:end); FIELD = FIELD(:,:,1:nskip:end);

switch ~FOURIER
    case 1 % Real space plot
        XNAME = '$x$'; YNAME = '$y$'; plane ='xy';
        [X,Y] = meshgrid(data.grids.x,data.grids.y);
        INTERP = 1;
        process = @(x) real(fftshift(ifourier_GENE(x)));
        shift_x = @(x) x;
        shift_y = @(x) x;
    case 0 % Frequencies plot
        XNAME = '$k_x$'; YNAME = '$k_y$'; plane ='kxky';
        [X,Y] = meshgrid(data.grids.kx,data.grids.ky);
        INTERP = 0;
        process = @(x) abs(fftshift(x,2));
        shift_x = @(x) fftshift(x,2);
        shift_y = @(x) fftshift(x,2); 
end

for it = 1:numel(TIME)
    tmp = squeeze(OPE_.*FIELD(:,:,it));
    F(:,:,it) = squeeze(process(tmp));
end

X         = shift_x(X);
Y         = shift_y(Y);
FILENAME  = [vname,'_',plane,'_obmp'];
FILENAME  = [DATADIR,FILENAME,format];
FIELDNAME = ['$',ltxname,'$'];
MOVIE.F   = F;
MOVIE.X   = X;
MOVIE.Y   = Y;
MOVIE.T   = TIME;
MOVIE.XNAME = XNAME;
MOVIE.YNAME = YNAME;
MOVIE.FPS = 16;
MOVIE.BWR = 0;
MOVIE.SAT = 0.75;
MOVIE.FIELDNAME = FIELDNAME;
MOVIE.FILENAME = FILENAME;
% Shoot the movie
if PLAY
    Ursula_Meier(MOVIE)
end

end