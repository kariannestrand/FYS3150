import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa

M = 7
T = 0.002
dt = 2.5e-5
N = int(T/dt)

probability_deviation = False
probability_time_evolution = True

if probability_deviation:
    U_cube_0 = pa.cx_cube()
    U_cube_0.load("v0_0.bin")

    U_cube_10 = pa.cx_cube()
    U_cube_10.load("v0_10.bin")

    t = np.linspace(0, T, N)

    p_0_array = np.zeros(N)
    p_10_array = np.zeros(N)
    for n in range(N):
        p_0 = 0.0
        p_10 = 0.0
        for i in range(M-2):
            for j in range(M-2):
                p_0 += (np.conj(U_cube_0[i, j, n])*U_cube_0[i, j, n]).real
                p_10 += (np.conj(U_cube_10[i, j, n])*U_cube_10[i, j, n]).real
        p_0_array[n] = p_0
        p_10_array[n] = p_10

    plt.plot(t, abs(1 - p_0_array), ".", label = "$v_0 = 0$")
    plt.plot(t, abs(1 - p_10_array), ".", label = "$v_0 = 10^{10}$")
    plt.title("Deviation from the normalized probability as a function of time", size = 22)
    plt.xlabel("Time", size = 22)
    plt.ylabel("|1 - $\sum_{i, j}\ p_{ij}$|", size = 22)
    plt.xticks(size = 20)
    plt.yticks(size = 20); plt.yscale("log")
    plt.legend(fontsize = 15)
    plt.show()


if probability_time_evolution:
    U_cube = pa.cx_cube()
    U_cube.load("T_0002.bin")

    U_0_Re = np.zeros((M-2, M-2))
    U_0001_Re = np.zeros((M-2, M-2))
    U_0002_Re = np.zeros((M-2, M-2))

    U_0_Im = np.zeros((M-2, M-2))
    U_0001_Im = np.zeros((M-2, M-2))
    U_0002_Im = np.zeros((M-2, M-2))

    p_0 = np.zeros((M-2, M-2))
    p_0001 = np.zeros((M-2, M-2))
    p_0002 = np.zeros((M-2, M-2))


    for i in range(M-2):
        for j in range(M-2):
            U_0_Re[i, j] = U_cube[i, j, 0].real
            U_0001_Re[i, j] = U_cube[i, j, int(N/2)].real
            U_0002_Re[i, j] = U_cube[i, j, N-1].real

            U_0_Im[i, j] = U_cube[i, j, 0].imag
            U_0001_Im[i, j] = U_cube[i, j, int(N/2)].imag
            U_0002_Im[i, j] = U_cube[i, j, N-1].imag

            p_0[i, j] = (np.conj(U_cube[i, j, 0])*U_cube[i, j, 0]).real
            p_0001[i, j] = (np.conj(U_cube[i, j, int(N/2)])*U_cube[i, j, int(N/2)]).real
            p_0002[i, j] = (np.conj(U_cube[i, j, N-1])*U_cube[i, j, N-1]).real



    fig = plt.figure(figsize = (8,6))

    p = [p_0, p_0001, p_0002]
    p_title = ["$p_0$", "$p_{0001}$", "$p_{0002}$"]

    U_Re = [U_0_Re, U_0001_Re, U_0002_Re]
    U_Re_title = ["Re($U_0$)", "Re($U_{0001}$)", "Re($U_{0002}$)"]

    U_Im = [U_0_Im, U_0001_Im, U_0002_Im]
    U_Im_title = ["Im($U_0$)", "Im($U_{0001}$)", "Im($U_{0002}$)"]

    for i in range(len(p)):
        plt.imshow(p[i])
        plt.title(p_title[i])
        plt.show()

    for i in range(len(p)):
        plt.imshow(U_Re[i])
        plt.title(U_Re_title[i])
        plt.show()

    for i in range(len(p)):
        plt.imshow(U_Im[i])
        plt.title(U_Im_title[i])
        plt.show()

"""
#
# Let's generate a dummy time series for a function z(x,y,t)
#

# Set up a 2D xy grid
h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)

# Array of time points
dt = 0.005
t_points = np.arange(0, 1+dt, dt)

# A function for a Gaussian that is travelling
# in the x direction and broadening as time passes
def z(x,y,t):
    v = 0.5
    x_c = 0.2
    sigma_x = 0.025 + 0.15 * t
    return 1. / (2 * np.pi * np.sqrt(sigma_x)) * np.exp(-0.5 * (x - x_c - v * t)**2 / sigma_x**2)

# Fill z_data_list with f(x,y,t)
z_data_list = []
for t in t_points:
    z_data = z(x, y, t)
    z_data_list.append(z_data)


#
# Now the list z_data_list contains a series of "frames" of z(x,y,t),
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

# Some settings
fontsize = 12
t_min = t_points[0]
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

# Plot the first frame
img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("z(x,y,t)", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(z_data_list[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)

# Run the animation!
plt.show()

# # Save the animation
# anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
"""
