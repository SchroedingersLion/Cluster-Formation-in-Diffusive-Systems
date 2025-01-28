import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import numpy as np
import os
import argparse
os.chdir(os.path.abspath(os.path.dirname(__file__)))

mpl.rcParams["savefig.directory"] = ""
plt.rcParams.update({'font.size': 30})
plt.rc('legend', fontsize=25)
plt.rc('figure', titlesize=36)
plt.rc('axes', titlesize=30) 
plt.rcParams['axes.grid'] = True
plt.rc('lines', linewidth=4)
plt.rcParams['axes.titlepad'] = 26
plt.rcParams['axes.labelsize'] = 30
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22) 

### Function to update animation.
def update(frame, x_coords, y_coords, z_coords, times, scatterplot, time_text):
    
    # for each frame, update the data stored on each artist.
    x = x_coords[frame]
    y = y_coords[frame]
    z = z_coords[frame]

    scatterplot._offsets3d=(x,y,z)


    # for i in range(0,len(x)):
    #     xdata = x[i]
    #     ydata = y[i]
    #     zdata = z[i]
    #     data = np.stack([xdata, ydata, zdata]).T
    #     scatterplot._offsets3d(data)
    #     scatterplot.set_offsets(data)

    time_text.set_text(f'Time: {times[frame]}')

    return (scatterplot)


# Parse the command-line arguments.
parser = argparse.ArgumentParser(description='Reads a trajectory file obtained by IPS.exe and creates an animation.')

parser.add_argument('file', type=str, help='Trajectory data file to animate.')
parser.add_argument('--interval', type=int, default=10, help='Time in ms between two frames in animation.')
parser.add_argument('--title', type=str, default="", help='Title of the animation.')
parser.add_argument('--save', type=str, default="", help='File name for saving the animation. Default: File not saved.')
parser.add_argument('--fps', type=int, default=3, help='Frames per second in the saved animation (overrides "interval").')

args = parser.parse_args()

filename = args.file
interval = args.interval
title = args.title
save_file = args.save
fps = args.fps

# filename = "3D_simu.csv_trajectory"
# interval = 10
# title = "test 3D simulation"
# save_file = "save_test.gif"
# fps = 3


# Load trajectory data.
arr = np.loadtxt(filename, delimiter=" ", skiprows=1) 

# Determine boxsize from coordinates.
xlim = [np.min(arr[:,1]), np.max(arr[:,1])]
ylim = [np.min(arr[:,2]), np.max(arr[:,2])]
zlim = [np.min(arr[:,3]), np.max(arr[:,3])]

# Create animation window.
# fig, ax = plt.subplots(figsize=(10, 10))
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set(xlim=xlim, ylim=ylim, zlim=zlim)
ax.set_xlabel('x', labelpad=15)
ax.set_ylabel('y', labelpad=15)
ax.set_zlabel('z', labelpad=15)
ax.set_title(title)

# Configure 3D plot
ax.set_proj_type('persp')  # Set perspective projection
ax.grid(True)  # Remove 2D grid
ax.set_box_aspect([1, 1, 1])  # Ensure the 3D aspect ratio is correct

# Compute number of particles (given by the first row where the time is no longer 0).
N_part = np.where(arr[:,0]!=0)[0][0]

# Create array of frames.
N_frames = len(arr)//N_part                      # Number of frames.
frames = arr.reshape(N_frames, N_part, 4)       

# Extract coordinates to feed into animation function.
x_coords = np.zeros((N_frames, N_part))
y_coords = np.zeros((N_frames, N_part))
z_coords = np.zeros((N_frames, N_part))
times = np.zeros((N_frames))
for i in range(0, N_frames):
    x_coords[i,:] = arr[i*N_part : (i+1)*N_part, 1]
    y_coords[i,:] = arr[i*N_part : (i+1)*N_part, 2]
    z_coords[i,:] = arr[i*N_part : (i+1)*N_part, 3]
    times[i] = arr[i*N_part, 0]

# Create scatter plot for animation function.
scat = ax.scatter(x_coords[0], y_coords[0], z_coords[0], c="b", s=5,)
time_text = ax.text2D(0.05, 0.9, "", transform=ax.transAxes)

# Plot animation.
ani = animation.FuncAnimation(fig=fig, func=update, frames=N_frames, fargs=(x_coords, y_coords, z_coords, times, scat, time_text), interval=interval)
# Save animation.
if save_file:
    print("Saving animation...")
    ani.save(save_file, writer='ffmpeg', fps=fps)

plt.show()


