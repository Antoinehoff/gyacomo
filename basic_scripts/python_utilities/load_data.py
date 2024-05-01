import h5py
import numpy as np
import tools
def load_data_0D(filename,dname):
    with h5py.File(filename, 'r') as file:
        # Load time data
        time  = file['data/var0d/time']
        time  = time[:]
        var0D = file['data/var0d/'+dname]
        var0D = var0D[:]
    return time, var0D

def load_data_3D_frame(filename,dname,tframe):
    with h5py.File(filename, 'r') as file:
        # load time
        time  = file['data/var3d/time']
        time  = time[:]
        # find frame
        iframe = tools.closest_index(time,tframe)
        tf     = time[iframe]
        # Load data
        if dname == 'Ni00':
            data = file[f'data/var3d/Na00/{iframe:06d}']
            data = data[:,:,:,0]
        else:
            data = file[f'data/var3d/{dname}/{iframe:06d}']
        data = data[:]
        data = data['real']+1j*data['imaginary'] 
        # Load the grids
        kxgrid = file[f'data/grid/coordkx']
        kygrid = file[f'data/grid/coordky']
        zgrid  = file[f'data/grid/coordz']
    return time,data,tf

def load_grids(filename):
    with h5py.File(filename, 'r') as file:
        # Load the grids
        kxgrid = file[f'data/grid/coordkx']
        kygrid = file[f'data/grid/coordky']
        zgrid  = file[f'data/grid/coordz']
        pgrid  = file[f'data/grid/coordp']
        jgrid  = file[f'data/grid/coordj']
        Nx     = kxgrid.size
        Nky    = kygrid.size
        Ny     = 2*(Nky-1)
        Lx     = 2*np.pi/kxgrid[1]
        Ly     = 2*np.pi/kygrid[1]
        xgrid  = np.linspace(-Lx/2,Lx/2,Nx)
        ygrid  = np.linspace(-Ly/2,Ly/2,Ny)
    return xgrid,kxgrid,ygrid,kygrid,zgrid,pgrid,jgrid

def load_group(filename,group):
    data = {}
    with h5py.File(filename, 'r') as file:
        g_  = file[f"data/"+group]
        for key in g_.keys():
            name='data/'+group+'/'+key
            data[key] = file.get(name)[:]
    return data

def load_params(filename):
    with h5py.File(filename, 'r') as file:
        nml_str = file[f"files/STDIN.00"][0]
        nml_str = nml_str.decode('utf-8')
        params = read_namelist(nml_str)
    return params

# Function to read all namelists from a file
def read_namelist(nml_str):
    Nspecies = 1
    all_namelists = {}
    current_namelist = None
    nml_str = nml_str.split('\n')
    for line in nml_str:
        line = line.split('!')
        line = line[0]
        if line.startswith('&'):
            current_namelist = line[1:].strip()
            if current_namelist == 'SPECIES':
                current_namelist = current_namelist + "_" + str(Nspecies)
                Nspecies = Nspecies + 1
            all_namelists[current_namelist] = {}
        elif line.startswith('/'):
            current_namelist = None
        elif current_namelist:
            parts = line.split('=')
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip().rstrip(',').strip("'").strip()
                if tools.is_convertible_to_float(value):
                    all_namelists[current_namelist][key] = float(value)
                else:
                    all_namelists[current_namelist][key] = value
    return all_namelists