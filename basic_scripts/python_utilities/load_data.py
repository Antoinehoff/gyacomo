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