from matplotlib.colors import Normalize
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa

h = 0.005
M = int(1.0/h + 1.0);
T = 0.01
dt = 2.5e-5
N = int(T/dt)

save_fig = True

probability_deviation = False
probability_deviation_bin_name_no_barrier = "probability_deviation_no_barrier.bin"
probability_deviation_bin_name_double_slit = "probability_deviation_double_slit.bin"

probability_time_evolution = False
probability_time_evolution_bin_name = "probability_time_evolution.bin"

detection_probability = False
detection_probability_bin_name_single = "detection_probability_single.bin"
detection_probability_bin_name_double = "detection_probability_double.bin"
detection_probability_bin_name_triple = "detection_probability_triple.bin"

animation = True
animation_bin_name = "animation.bin"


if probability_deviation:
    U_cube_0_T = pa.cx_cube()
    U_cube_0_T.load(probability_deviation_bin_name_no_barrier)
    U_cube_0 = np.transpose(U_cube_0_T)


    U_cube_10_T = pa.cx_cube()
    U_cube_10_T.load(probability_deviation_bin_name_double_slit)
    U_cube_10 = np.transpose(U_cube_10_T)


    p_0_array = np.zeros(N)
    p_10_array = np.zeros(N)
    for n in range(N):
        p_0 = 0.0
        p_10 = 0.0
        for i in range(M-2):
            for j in range(M-2):
                p_0 += np.real(np.conj(U_cube_0[i, j, n])*U_cube_0[i, j, n])
                p_10 += np.real(np.conj(U_cube_10[i, j, n])*U_cube_10[i, j, n])
        p_0_array[n] = p_0
        p_10_array[n] = p_10


    t = np.linspace(0, T, N)
    plt.figure(figsize=(11,8))
    plt.plot(t, abs(1 - p_0_array), ".", label = "$v_0 = 0$")
    plt.plot(t, abs(1 - p_10_array), ".", label = "$v_0 = 10^{10}$")
    plt.xlabel("Time", size = 24)
    plt.ylabel("|1 - $\sum_{i, j}\ p_{ij}$|", size = 22)
    plt.xticks(size = 24)
    plt.yticks(size = 24); plt.yscale("log")
    plt.legend(fontsize = 24)
    plt.tight_layout()
    if save_fig:
        plt.savefig('probability_deviation.pdf')
    plt.show()

if probability_time_evolution:
    U_cube_T = pa.cx_cube()
    U_cube_T.load(probability_time_evolution_bin_name)
    U_cube = np.transpose(U_cube_T)

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

            p_0[i, j] = np.real(np.conj(U_cube[i, j, 0])*U_cube[i, j, 0])
            p_0001[i, j] = (np.conj(U_cube[i, j, int(N/2)])*U_cube[i, j, int(N/2)]).real
            p_0002[i, j] = (np.conj(U_cube[i, j, N-1])*U_cube[i, j, N-1]).real


    p = [p_0, p_0001, p_0002]
    p_title = ["$p_0$", "$p_{0001}$", "$p_{0002}$"]

    U_Re = [U_0_Re, U_0001_Re, U_0002_Re]
    U_Re_title = ["Re($U_0$)", "Re($U_{0001}$)", "Re($U_{0002}$)"]

    U_Im = [U_0_Im, U_0001_Im, U_0002_Im]
    U_Im_title = ["Im($U_0$)", "Im($U_{0001}$)", "Im($U_{0002}$)"]


    x = np.linspace(0+h, 1-h, M-2);
    y = np.linspace(0+h, 1-h, M-2);
    X, Y = np.meshgrid(x,y)

    t_string = ["0", "0.001", "0.002"]

    for i in range(len(p)):
        plt.contourf(X, Y, np.sqrt(p[i]), 20)
        plt.xlabel("x", size = 17)
        plt.ylabel("y", size = 17)
        plt.xticks(size = 15)
        plt.yticks(size = 15)
        cb = plt.colorbar()
        cb.set_label(label = 'p(x, y; t = ' + t_string[i] + '$)^{1/2}$', size = 17)
        cb.ax.tick_params(labelsize = 15)
        if save_fig:
            plt.savefig('prob' + t_string[i] + '.pdf')
        plt.show()
        plt.close("all")


    for i in range(len(p)):
        plt.contourf(X, Y, U_Re[i], 20)
        plt.xlabel("x", size = 17)
        plt.ylabel("y", size = 17)
        plt.xticks(size = 15)
        plt.yticks(size = 15)
        cb = plt.colorbar()
        cb.set_label(label='Re(u(x, y, t = ' + t_string[i] + '))', size = 17)
        cb.ax.tick_params(labelsize = 15)
        if save_fig:
            plt.savefig('real' + t_string[i] + '.pdf')
        plt.show()
        plt.close("all")

    for i in range(len(p)):
        plt.contourf(X, Y, U_Im[i], 20)
        plt.xlabel("x", size = 17)
        plt.ylabel("y", size = 17)
        plt.xticks(size = 15)
        plt.yticks(size = 15)
        cb = plt.colorbar()
        cb.set_label(label='Im(u(x, y, t = ' + t_string[i] + '))', size = 17)
        cb.ax.tick_params(labelsize = 14)
        if save_fig:
            plt.savefig('imag' + t_string[i] + '.pdf')
        plt.show()
        plt.close("all")

if detection_probability:
    U_cube_single_T = pa.cx_cube()
    U_cube_single_T.load(detection_probability_bin_name_single)
    U_cube_single = np.transpose(U_cube_single_T)

    U_cube_double_T = pa.cx_cube()
    U_cube_double_T.load(detection_probability_bin_name_double)
    U_cube_double = np.transpose(U_cube_double_T)

    U_cube_triple_T = pa.cx_cube()
    U_cube_triple_T.load(detection_probability_bin_name_triple)
    U_cube_triple = np.transpose(U_cube_triple_T)

    p_single = np.zeros((M-2, M-2))
    p_double = np.zeros((M-2, M-2))
    p_triple = np.zeros((M-2, M-2))

    for i in range(M-2):
        for j in range(M-2):
            p_single[i, j] = (np.conj(U_cube_single[i, j, N-1])*U_cube_single[i, j, N-1]).real
            p_double[i, j] = (np.conj(U_cube_double[i, j, N-1])*U_cube_double[i, j, N-1]).real
            p_triple[i, j] = (np.conj(U_cube_triple[i, j, N-1])*U_cube_triple[i, j, N-1]).real

    x = np.linspace(0+h, 1-h, M-2);
    y = np.linspace(0+h, 1-h, M-2);
    x_index = int(np.where(x==0.8)[0])

    plt.figure(figsize=(11,8))
    plt.plot(y, p_single[:, x_index], label = "Single-slit experiment")
    plt.plot(y, p_double[:, x_index], label = "Double-slit experiment")
    plt.plot(y, p_triple[:, x_index], label = "Triple-slit experiment")
    plt.xlabel("y", size = 24)
    plt.ylabel("p(y | x = 0.8; t = 0.002)", size = 22)
    plt.xticks(size = 24)
    plt.yticks(size = 24);
    plt.legend(fontsize = 17)
    plt.tight_layout()
    if save_fig:
        plt.savefig('detection_probability.pdf')
    plt.show()

if animation:
    U_cube_T = pa.cx_cube()
    U_cube_T.load(animation_bin_name)
    U_cube = np.transpose(U_cube_T)

    x_points = np.arange(0, 1+h, h);
    y_points = np.arange(0, 1+h, h);
    x, y = np.meshgrid(x_points, y_points, sparse = True)

    t_points = np.arange(0, 1+dt, dt)

    p_data_list = []
    for n in range(N):
        p_data = np.real(np.conj(U_cube[:, :, n])*U_cube[:, :, n])
        p_data_list.append(p_data)


    fontsize = 14
    t_min = t_points[0]
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(p_data_list[0]))

    # Plot the first frame
    img = ax.imshow(p_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("p(x, y, t)", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(p_data_list[i]))
        img.set_norm(norm)

        # Update z data
        img.set_data(p_data_list[i])

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(p_data_list), 2), repeat=False, blit=0)

    # Run the animation!
    plt.show()

    # # Save the animation
    if save_fig:
        anim.save('animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
