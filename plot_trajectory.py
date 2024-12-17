import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
#plt.rc('ytick', labelsize=22) 




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


# Load trajectory data.
arr = np.loadtxt(filename, delimiter=" ", skiprows=1) 

# Determine boxsize from coordinates.
ylim = [np.min(arr[:,1]), np.max(arr[:,1])]
threshold = (ylim[1]-ylim[0])/2 # Half of the boxsize, used to exclude periodic boundary jumps from plot.

# Create animation window.
fig, ax = plt.subplots(figsize=(10, 10))
ax.set(ylim=ylim, xlabel='Time', ylabel='x')
ax.set_title(title)
colors = mpl.cm.tab10.colors

# Compute number of particles (given by the first row where the time is no longer 0).
N_part = np.where(arr[:,0]!=0)[0][0]

# Create array of frames.
N_frames = len(arr)//N_part                      # Number of frames.
frames = arr.reshape(N_frames, N_part, 2)       

# Extract coordinates to feed into animation function.
x_coords = np.zeros((N_frames, N_part))

times = np.zeros((N_frames))
for i in range(0, N_frames):
    x_coords[i,:] = arr[i*N_part : (i+1)*N_part, 1]
    times[i] = arr[i*N_part, 0]

# Plot particle trajectories one by one.
for i in range(0, N_part):

    positions = x_coords[:,i]

    # Need to split the data to avoid plotting of jumps caused by periodic boundaries.
    jumps = np.abs(np.diff(positions)) > threshold
    # Create segments for continuous plotting
    segments = []
    start_idx = 0
    for k, jump in enumerate(jumps):
        if jump:  # End current segment at the jump
            segments.append((times[start_idx:k + 1], positions[start_idx:k + 1]))
            start_idx = k + 1
    # Append the last segment
    segments.append((times[start_idx:], positions[start_idx:]))
    
    color = colors[i % len(colors)]
    for t_seg, p_seg in segments:
        ax.plot(t_seg, p_seg, color=color)

# Save animation.
# if save_file:
#     print("Saving animation...")
#     ani.save(save_file, writer='ffmpeg', fps=fps)

plt.show()


