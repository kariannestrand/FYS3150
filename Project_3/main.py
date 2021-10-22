import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pyarma as pa


n = 1                            # number of particles
interaction = False

z_t = False
x_y = False
phase_space = False
trajectory = False
relative_error = False
error_convergence_rate = False

v_x = phase_space
v_y = phase_space
v_z = phase_space


# constants
q = 1.
m = 40.078

B0 = 9.65e1
V0 = 9.65e8
d = 1.0e4

omega_0 = q*B0/m
omega_z = np.sqrt(2*q*V0/(m*d**2))

omega_p = 0.5*(omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2));
omega_m = 0.5*(omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2));


def V(filename_v):
    V = pa.mat()
    V.load(filename_v)
    vel = np.transpose(V)

    x = vel[:, 0]
    y = vel[:, 1]
    z = vel[:, 2]

    return x, y, z


def R(filename_r):
    R = pa.mat()
    R.load(filename_r)
    pos = np.transpose(R)

    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]

    return x, y, z


def r_analytical(filename_r, filename_v, h):
    Rx = R(filename_r)[0]
    Rz = R(filename_r)[2]

    Vy = V(filename_v)[1]

    N = len(Rx)

    x_0 = Rx[0]
    v_0 = Vy[0]
    z_0 = Rz[0]


    Ap = (v_0 + omega_m*x_0)/(omega_m - omega_p)
    Am = - (v_0 + omega_p*x_0)/(omega_m - omega_p)


    # analytical solution, r_exact
    x = np.empty(N)
    y = np.empty(N)
    z = np.empty(N)

    r_exact = np.empty(N)

    for i in range(N):
        x[i] = Ap*np.cos(omega_p*i*h) + Am*np.cos(omega_m*i*h)
        y[i] = - Ap*np.sin(omega_p*i*h) - Am*np.sin(omega_m*i*h)
        z[i] = z_0*np.cos(omega_z*i*h)

        r_exact[i] = np.sqrt(x[i]**2 + y[i]**2 + z[i]**2)

    return r_exact


def r_numerical(filename_r):
    Rx = R(filename_r)[0]
    Ry = R(filename_r)[1]
    Rz = R(filename_r)[2]

    N = len(Rx)

    r_num = np.empty(N)
    for i in range(N):
        r_num[i] = np.sqrt(Rx[i]**2 + Ry[i]**2 + Rz[i]**2)

    return r_num


def relative_error(filename_r, filename_v, h):
    r_exact = r_analytical(filename_r, filename_v, h)
    r_num = r_numerical(filename_r)

    rel_err = np.abs((r_exact - r_num)/r_exact)

    return rel_err

def error_convergence_rate(filename_r_list, filename_v_list, h_list):
    N = len(h_list)
    delta_max = np.empty(N)

    for i in range(N):
        r_exact = r_analytical(filename_r_list[i], filename_v_list[i], h_list[i])
        r_num = r_numerical(filename_r_list[i])

        delta_max[i] = np.max(np.abs(r_exact - r_num))

    r_err = 0
    for i in range(1, N):
        r_err += 0.25*np.log10(delta_max[i]/delta_max[i-1])/np.log10(h_list[i]/h_list[i-1])

    return r_err

def time(filename_r):
    N = len(r_numerical(filename_r))
    t = np.linspace(0, 100, N)
    return t


filename_r = "r_0_0001.bin"
filename_v = "v_0_0001.bin"

if z_t:
    #fig = plt.figure()
    Rz = R(filename_r)[2]
    t = time(filename_r)

    plt.plot(t, Rz)

    if n == 1:
        plt.title("Movement of one particle in the z-direction", size = 12)

    if n == 2:
        plt.plot(t, Rz2)

        if interaction:
            plt.title("Movement of two particles in the z-direction w/ interaction", size = 12)
        else:
            plt.title("Movement of two particles in the z-direction w/o interaction", size = 12)

    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("z/[$\mu$m]", size = 12)

    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    if n == 1:
        plt.savefig('pdf/1zt.pdf')
    if n == 2:
        if interaction:
            plt.savefig('pdf/zt_2_w.pdf')
        else:
            plt.savefig('pdf/zt_2_wo.pdf')

    plt.show()

if x_y:
    #fig = plt.figure()
    Rx = R(filename_r)[0]
    Ry = R(filename_r)[1]
    plt.plot(Rx, Ry)

    if n == 1:
        plt.title("Movement of one particle in the xy-plane", size = 10)

    if n == 2:
        plt.plot(Rx2, Ry2)

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
            plt.savefig('pdf/xy_2_w.pdf')
        else:
            plt.savefig('pdf/xy_2_wo.pdf')


    plt.show()

if v_x:
    #fig = plt.figure()
    Rx = R(filename_r)[0]
    Vx = V(filename_v)[0]
    plt.plot(Rx, Vx)

    if n == 1:
        plt.title("Velocity as function of position in the x-direction", size = 10)

    if n == 2:
        plt.plot(Rx2, Vx2)

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
    #fig = plt.figure()
    Ry = R(filename_r)[1]
    Vy = V(filename_v)[1]
    plt.plot(Ry, Vy)

    if n == 1:
        plt.title("Velocity as function of position in the y-direction", size = 10)

    if n == 2:
        plt.plot(Ry2, Vy2)

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
    #fig = plt.figure()
    Rz = R(filename_r)[2]
    Vz = V(filename_v)[2]
    plt.plot(Rz, Vz)

    if n == 1:
        plt.title("Velocity as function of position in the z-direction", size = 10)

    if n == 2:
        plt.plot(Rz2, Vz2)

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
            plt.savefig('pdf/vz_2_w.pdf')
        else:
            plt.savefig('pdf/vz_2_wo.pdf')

    plt.show()

if trajectory:
    #fig = plt.figure()

    Rx = R(filename_r)[0]
    Ry = R(filename_r)[1]
    Rz = R(filename_r)[2]

    ax = plt.axes(projection="3d")
    plt.tight_layout()
    ax.plot3D(Rx, Ry, Rz, 'blue', label='Trajectory of particle 1')
    plt.plot(Rx[0], Ry[0], Rz[0], "ro")

    if n == 1:
        plt.title("Trajectory of one particle", size = 12)

    if n == 2:
        ax.plot3D(Rx2, Ry2, Rz2, 'red', label ='Trajectory of particle 2')
        plt.plot(Rx2[0], Ry2[0], Rz2[0], "ro")

        if interaction:
            plt.title("Trajectory of two particles w/ interaction", size = 12)
        else:
            plt.title("Trajectory of two particles w/o interaction", size = 12)

    ax.set_xlabel("x/[$\mu$m]", size = 12)
    ax.set_ylabel("y/[$\mu$m]", size = 12)
    ax.set_zlabel("z/[$\mu$m]", size = 12)

    plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

    plt.legend()
    if n == 1:
        plt.savefig('pdf/3D_1.pdf')

    if n == 2:
        if interaction:
            plt.savefig('pdf/3D_2_w.pdf')
        else:
            plt.savefig('pdf/3D_2_wo.pdf')

    plt.show()


if relative_error:
    filename_r_list = ["r_0_0001.bin", "r_0_0001.bin", "r_0_0001.bin", "r_0_0001.bin", "r_0_0001.bin"]
    filename_v_list = ["v_0_0001.bin", "v_0_0001.bin", "v_0_0001.bin", "v_0_0001.bin", "v_0_0001.bin"]
    h_list = [1., 0.1, 0.01, 0.001, 0.0001]


    N = len(h_list)
    for i in range(N):
        t = time(filename_r_list[i])
        rel_err = relative_error(filename_r_list[i], filename_v_list[i], h_list[i])
        plt.plot(t, rel_err, label = "h = " + str(h_list[i]))

    plt.title("Relative Error with RK4", size = 12)
    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("$|(r_{exact} - r_{numerical})/r_{numerical}|$", size = 12)
    plt.legend()
    plt.savefig('pdf/rel_err_RK4.pdf')
    plt.show()

if error_convergence_rate:
    filename_r_list = ["r_0_0001.bin", "r_0_0001.bin", "r_0_0001.bin", "r_0_0001.bin", "r_0_0001.bin"]
    filename_v_list = ["v_0_0001.bin", "v_0_0001.bin", "v_0_0001.bin", "v_0_0001.bin", "v_0_0001.bin"]
    h_list = [1., 0.1, 0.01, 0.001, 0.0001]

    r_err = error_convergence_rate(filename_r_list, filename_v_list, h_list)
    print("Error convergence rate with RK4: r_err = {}".format(r_err))
