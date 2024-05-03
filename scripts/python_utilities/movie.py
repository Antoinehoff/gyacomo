import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import load_data as loader
from matplotlib.widgets import Slider
import sys
import tools

def update_mov(frame):
    # plt.imshow(data[:, :, frame], cmap='gray',interpolation='quadric')
    return [data[:, :, frame]]
    
def update_slide(val):
    frame = int(slider.val)
    im.set_data(data[:, :, frame])
    fig.canvas.draw_idle()

# Parameters (default values)
option    = 'movie'
filename  = 'outputs_00.h5'
fieldname = 'Tpar' 
if len(sys.argv) > 1:
    option   = sys.argv[1]
    if option == 'help':
        print("This script produces a movie of a field in a 2D perpendicular plane")
        print("Usage: python3 movie.py [movie/slider (opt), outfilename (opt), phi/dens/Tpar etc. (opt)]")
        sys.exit(1)
if len(sys.argv) > 2:
    filename = sys.argv[2]
if len(sys.argv) > 3:
    fieldname = sys.argv[3]
if len(sys.argv) > 4:
    print("Usage: python3 movie.py [movie/slider (opt), outputs_XX.h5 (opt), phi/dens/Tpar etc. (opt)]")
    sys.exit(1)

# Load grids
x,kx,y,ky,z,p,j  = loader.load_grids(filename)
x0 = x[0]; x1 = x[-1]; y0 = y[0]; y1 = y[-1]
# Load first frame
t3D, field, tf  = loader.load_data_3D_frame(filename,fieldname,0)
field           = tools.zkxky_to_xy_const_z(field,-1)
shape = np.shape(field)
Nt    = np.size(t3D)
# Build the movie
data  = np.zeros([shape[0],shape[1],Nt])
it = 0
for t_ in t3D:
    t3D, field, tf  = loader.load_data_3D_frame(filename,fieldname,t_)
    field           = tools.zkxky_to_xy_const_z(field,-1)
    fmax            = np.max(np.abs(field))
    if fmax == 0: fmax = 1
    data[:,:,it]    = field/fmax
    it = it + 1

# Create figure and axis
fig, ax = plt.subplots()
plt.title(fieldname)
plt.xlabel(r'$x/\rho_s$')
plt.ylabel(r'$y/\rho_s$')

# Display the first frame
im = ax.imshow(data[:, :, 0], cmap='seismic',interpolation='quadric',extent=[x0,x1,y0,y1])
im.set_clim(-1, 1)

if option == 'movie' or option == 'save':
    # Define the animation
    # ani = animation.FuncAnimation(fig, update_mov, frames=data.shape[2], interval=50)
    ims = []
    for i in range(Nt):
        im = ax.imshow(data[:, :, i], cmap='seismic',interpolation='quadric',extent=[x0,x1,y0,y1])
        ims.append([im])
    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True)
elif option == 'slider':
    # Add slider widget
    slider_ax = plt.axes([0.1, 0.01, 0.8, 0.03])
    slider = Slider(slider_ax, 'Time', 0, data.shape[2]-1, valinit=0, valstep=t3D[1]-t3D[0])
    slider.on_changed(update_slide)
else:
    print("Usage: python3 movie.py [movie, slider (opt)]")
    sys.exit(1)

if option == 'save':
    ani.save('movie_'+fieldname+'.mov', fps=30)

# Show the animation
plt.show()