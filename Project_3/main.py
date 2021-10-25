import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pyarma as pa


n = 2                           # number of particles

interaction = False

z_t = False
x_y = False
phase_space = False
trajectory = True

particles_trapped = False

v_x = phase_space
v_y = phase_space
v_z = phase_space

def V(filename_v):
    V = pa.mat()
    V.load(filename_v)
    vel = np.transpose(V)

    x = vel[:, 0]
    y = vel[:, 1]
    z = vel[:, 2]

    return x, y, z

filename_v = 'RK4_v_0_interaction.bin'

def R(filename_r):
    R = pa.mat()
    R.load(filename_r)
    pos = np.transpose(R)

    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]

    return x, y, z


def time(filename_r):
    N = len(R(filename_r))
    t = np.linspace(0, 100, N)
    return t


def particles_trapped(filename):
    pt = np.loadtxt(filename, skiprows = 1)
    n_trapped = pt[:, 0]
    omega_v = pt[:, 1]

    return omega_v, n_trapped


filename_r = "RK4_r_0_0001dt.bin"
filename_v = "RK4_v_0_0001dt.bin"

if z_t:
    if n == 1:
        t = time(filename_r)
        Rz = R(filename_r)[2]
        plt.plot(t, Rz)
        plt.title("Movement of one particle in the z-direction", size = 12)

    if n == 2:
        filename_r_list = ["RK4_r_0_0001dt.bin", "RK4_r_1_0001dt.bin"]
        for i in range(n):
            t = time(filename_r_list[i])
            Rx = R(filename_r_list[i])[0]
            Ry = R(filename_r_list[i])[1]
            plt.plot(Rx, Ry)

            if interaction:
                plt.title("Movement of two particles in the z-direction w/ interaction", size = 12)
            else:
                plt.title("Movement of two particles in the z-direction w/o interaction", size = 12)


    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("z/[$\mu$m]", size = 12)

    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    if n == 1:
        plt.savefig('pdf/zt_1.pdf')
    if n == 2:
        if interaction:
            plt.savefig('pdf/zt_2_w.pdf')
        else:
            plt.savefig('pdf/zt_2_wo.pdf')

    plt.show()

if x_y:
    if n == 1:
        Rx = R(filename_r)[0]
        Ry = R(filename_r)[1]
        plt.plot(Rx, Ry)
        plt.title("Movement of one particle in the xy-plane", size = 10)

    if n == 2:
        filename_r_list = ["RK4_r_0_0001dt.bin", "RK4_r_1_0001dt.bin"]
        for i in range(n):
            Rx = R(filename_r_list[i])[0]
            Ry = R(filename_r_list[i])[1]
            plt.plot(Rx, Ry)

        if interaction:
            plt.title("Movement of two particles in the xy-plane w/ interaction", size = 10)
        else:
            plt.title("Movement of two particles in the xy-plane w/o interaction", size = 10)

    plt.xlabel("x/[$\mu$m]", size = 12)
    plt.ylabel("y/[$\mu$m]", size = 12)

    plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

    if n == 1:
        plt.savefig('pdf/xy_1.pdf')
    if n == 2:
        if interaction:
            plt.savefig('xy_2_w.pdf')
        else:
            plt.savefig('xy_2_wo.pdf')

    plt.show()


if v_x:
    if n == 1:
        Rx = R(filename_r)[0]
        Vx = V(filename_v)[0]
        plt.plot(Rx, Vx)
        plt.title("IDK YET", size = 10)

    if n == 2:
        filename_r_list = ["RK4_r_0_0001dt.bin", "RK4_r_1_0001dt.bin"]
        filename_v_list = ["RK4_v_0_0001dt.bin", "RK4_v_1_0001dt.bin"]
        for i in range(n):
            Rx = R(filename_r_list[i])[0]
            Vx = V(filename_v_list[i])[0]
            plt.plot(Rx, Ry)

        if interaction:
            plt.title("IDK YET", size = 10)
        else:
            plt.title("IDK YET", size = 10)

    plt.xlabel("x/[$\mu$m]", size = 12)
    plt.ylabel("v$_x$/[$\mu$m/$\mu$s]", size = 12)

    plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

    if n == 1:
        plt.savefig('pdf/vx_1.pdf')
    if n == 2:
        if interaction:
            plt.savefig('pdf/vx_2_w.pdf')
        else:
            plt.savefig('pdf/vx_2_wo.pdf')

    plt.show()

if v_y:
    if n == 1:
        Ry = R(filename_r)[1]
        Vy = V(filename_v)[1]
        plt.plot(Ry, Vy)
        plt.title("IDK YET", size = 10)

    if n == 2:
        filename_r_list = ["RK4_r_0_0001dt.bin", "RK4_r_1_0001dt.bin"]
        filename_v_list = ["RK4_v_0_0001dt.bin", "RK4_v_1_0001dt.bin"]
        for i in range(n):
            Ry = R(filename_r_list[i])[1]
            Vy = V(filename_v_list[i])[1]
            plt.plot(Ry, Vy)

        if interaction:
            plt.title("IDK YET", size = 10)
        else:
            plt.title("IDK YET", size = 10)

    plt.xlabel("y/[$\mu$m]", size = 12)
    plt.ylabel("v$_y$/[$\mu$m/$\mu$s]", size = 12)

    plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

    if n == 1:
        plt.savefig('pdf/vy_1.pdf')
    if n == 2:
        if interaction:
            plt.savefig('pdf/vy_2_w.pdf')
        else:
            plt.savefig('pdf/vy_2_wo.pdf')
    plt.show()

if v_z:
    if n == 1:
        Rz = R(filename_r)[2]
        Vz = V(filename_v)[2]
        plt.plot(Rz, Vz)
        plt.title("IDK YET", size = 10)

    if n == 2:
        #filename_r_list = ["RK4_r_0_0001dt.bin", "RK4_r_1_0001dt.bin"]
        #filename_v_list = ["RK4_v_0_0001dt.bin", "RK4_v_1_0001dt.bin"]
        filename_r_list = ["RK4_r_0_interaction.bin", "RK4_r_1_interaction.bin"]
        filename_v_list = ["RK4_v_0_interaction.bin", "RK4_v_1_interaction.bin"]
        for i in range(n):
            Rz = R(filename_r_list[i])[2]
            Vz = V(filename_v_list[i])[2]
            plt.plot(Rz, Vz)

        if interaction:
            plt.title("IDK YET", size = 10)
        else:
            plt.title("IDK YET", size = 10)

    plt.xlabel("z/[$\mu$m]", size = 12)
    plt.ylabel("v$_z$/[$\mu$m/$\mu$s]", size = 12)

    plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

    if n == 1:
        plt.savefig('pdf/vz_1.pdf')
    if n == 2:
        if interaction:
            plt.savefig('vz_2_w.pdf')
        else:
            plt.savefig('vz_2_wo.pdf')

    plt.show()


if trajectory:
    ax = plt.axes(projection="3d")
    plt.tight_layout()

    if n == 1:
        Rx = R(filename_r)[0]
        Ry = R(filename_r)[1]
        Rz = R(filename_r)[2]

        ax.plot3D(Rx, Ry, Rz, 'blue', label='Trajectory of particle 1')
        plt.plot(Rx[0], Ry[0], Rz[0], "ro")


    if n == 2:
        filename_r_list = ["RK4_r_0_0001dt.bin", "RK4_r_1_0001dt.bin"]
        for i in range(n):
            Rx = R(filename_r_list[i])[0]
            Ry = R(filename_r_list[i])[1]
            Rz = R(filename_r_list[i])[2]

            ax.plot3D(Rx, Ry, Rz, label='Trajectory of particle ' + str(i+1))


    ax.set_xlabel("x/[$\mu$m]", size = 12)
    ax.set_ylabel("y/[$\mu$m]", size = 12)
    ax.set_zlabel("z/[$\mu$m]", size = 12)

    plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

    plt.legend()
    if n == 1:
        plt.savefig('pdf/3D_1.pdf')

    if n == 2:
        if interaction:
            plt.savefig('finished_plots/3D_2_w.pdf')
        else:
            plt.savefig('finished_plots/3D_2_wo.pdf')

    plt.show()


if particles_trapped:
    filename_list = ["trapped_01f.txt", "trapped_04f.txt", "trapped_07f.txt"]
    f_list = [0.1, 0.4, 0.7]

    for i in range(len(filename_list)):
        filename = filename_list[i]
        omega_v = particles_trapped(filename)[0]
        n_trapped = particles_trapped(filename)[1]
        plt.plot(omega_v, n_trapped/n, label = "f = " + str(f_list[i]))

    plt.legend()
    plt.show()


