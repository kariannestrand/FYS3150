import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


n = 1                            # number of particles

# buttons
interaction = False
read_velocity = True             # must always be True
read_position = True             # must always be True

z_t = False
x_y = False

phase_space = False
trajectory = False
relative_error = True
error_convergence_rate = True



v_x = phase_space
v_y = phase_space
v_z = phase_space

relative_error_RK4 = relative_error
relative_error_Euler = relative_error

error_convergence_rate_RK4 = error_convergence_rate
error_convergence_rate_Euler = error_convergence_rate



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


def r_analytical(filename_r, filename_v, h):

    pos = np.loadtxt(filename_r, skiprows = 1)
    Rx1 = pos[:, 0]
    Rz1 = pos[:, 2]

    vel = np.loadtxt(filename_v, skiprows = 1)
    Vy1 = pos[:, 1]

    N = len(Rx1)

    x_0 = Rx1[0]
    v_0 = Vy1[0]
    z_0 = Rz1[0]


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
    pos = np.loadtxt(filename_r, skiprows = 1)
    Rx1 = pos[:, 0]
    Ry1 = pos[:, 1]
    Rz1 = pos[:, 2]

    N = len(Rx1)

    r_num = np.empty(N)
    for i in range(N):
        r_num[i] = np.sqrt(Rx1[i]**2 + Ry1[i]**2 + Rz1[i]**2)

    """
    t = np.linspace(0, 100, N)
    plt.plot(t, r_num)
    plt.show()
    """

    return r_num


def relative_error(filename_r, filename_v, h):
    r_exact = r_analytical(filename_r, filename_v, h)
    r_num = r_numerical(filename_r)

    rel_err = np.abs((r_exact - r_num)/r_exact)

    return rel_err


def time(filename_r):
    N = len(r_numerical(filename_r))
    t = np.linspace(0, 100, N)
    return t


if read_velocity:
    filename = "RK4_v_1_0001dt.txt"
    vel = np.loadtxt(filename, skiprows = 1)
    Vx1 = vel[:, 0]
    Vy1 = vel[:, 1]
    Vz1 = vel[:, 2]


    if n == 2:
        filename = "RK4_v_2_0001dt.txt"
        vel = np.loadtxt(filename, skiprows = 1)
        Vx2 = vel[:, 0]
        Vy2 = vel[:, 1]
        Vz2 = vel[:, 2]


if read_position:
    filename = "RK4_r_1_0001dt.txt"
    pos = np.loadtxt(filename, skiprows = 1)
    Rx1 = pos[:, 0]
    Ry1 = pos[:, 1]
    Rz1 = pos[:, 2]


    if n == 2:
        filename = "RK4_r_2_0001dt.txt"
        pos = np.loadtxt(filename, skiprows = 1)
        Rx2 = pos[:, 0]
        Ry2 = pos[:, 1]
        Rz2 = pos[:, 2]


if z_t:
    t = np.linspace(0, 100, len(Rz1))

    z_min = np.min(Rz1)
    index = np.where(Rz1 == Rz1.min())
    f = 1/(np.max(t[index]) - np.min(t[index]))

    print(omega_z)    # frequency of particle about the z-axis
    print(f)          # trying to calculate frequency from plot


    fig = plt.figure()
    plt.plot(t, Rz1)


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
    fig = plt.figure()
    plt.plot(Rx1, Ry1)

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
    fig = plt.figure()
    plt.plot(Rx1, Vx1)

    if n == 1:
        plt.title("IDK YET", size = 10)

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
    fig = plt.figure()
    plt.plot(Ry1, Vy1)

    if n == 1:
        plt.title("IDK YET", size = 10)

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
    fig = plt.figure()
    plt.plot(Rz1, Vz1)

    if n == 1:
        plt.title("IDK YET", size = 10)

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
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    plt.tight_layout()
    ax.plot3D(Rx1, Ry1, Rz1, 'blue', label='Trajectory of particle 1')
    plt.plot(Rx1[0], Ry1[0], Rz1[0], "ro")

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


if relative_error_RK4:
    filename_r_RK4 = ["RK4_r_1_1dt.txt", "RK4_r_1_01dt.txt", "RK4_r_1_001dt.txt", "RK4_r_1_0001dt.txt", "RK4_r_1_00001dt.txt"]
    filename_v_RK4 = ["RK4_v_1_1dt.txt", "RK4_v_1_01dt.txt", "RK4_v_1_001dt.txt", "RK4_v_1_0001dt.txt", "RK4_v_1_00001dt.txt"]

    h = [1., 0.1, 0.01, 0.001, 0.0001]
    for i in range(len(filename_r_RK4)):
        plt.plot(time(filename_r_RK4[i]), relative_error(filename_r_RK4[i], filename_v_RK4[i], h[i]), label = "h = " + str(h[i]))

    plt.title("Relative Error with RK4", size = 12)

    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("$|(r_{exact} - r_{numerical})/r_{numerical}|$", size = 12)

    #plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
    plt.legend()

    plt.savefig('pdf/rel_err_RK4.pdf')

    plt.show()


if relative_error_Euler:
    filename_r_Euler = ["Euler_r_1_1dt.txt", "Euler_r_1_01dt.txt", "Euler_r_1_001dt.txt", "Euler_r_1_0001dt.txt", "Euler_r_1_00001dt.txt"]
    filename_v_Euler = ["Euler_v_1_1dt.txt", "Euler_v_1_01dt.txt", "Euler_v_1_001dt.txt", "Euler_v_1_0001dt.txt", "Euler_v_1_00001dt.txt"]

    h = [1., 0.1, 0.01, 0.001, 0.0001]
    for i in range(len(filename_r_Euler)):
        plt.plot(time(filename_r_Euler[i]), relative_error(filename_r_Euler[i], filename_v_Euler[i], h[i]), label = "h = " + str(h[i]))

    plt.title("Relative Error with Forward Euler", size = 12)

    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("$|(r_{exact} - r_{numerical})/r_{numerical}|$", size = 12)

    #plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
    plt.legend()

    plt.savefig('pdf/rel_err_Euler.pdf')

    plt.show()


if error_convergence_rate_RK4:
    filename_r_RK4 = ["RK4_r_1_1dt.txt", "RK4_r_1_01dt.txt", "RK4_r_1_001dt.txt", "RK4_r_1_0001dt.txt", "RK4_r_1_00001dt.txt"]
    filename_v_RK4 = ["RK4_v_1_1dt.txt", "RK4_v_1_01dt.txt", "RK4_v_1_001dt.txt", "RK4_v_1_0001dt.txt", "RK4_v_1_00001dt.txt"]

    h = [1., 0.1, 0.01, 0.001, 0.0001]

    Delta_max = np.empty(len(h))
    for i in range(len(h)):
        r_exact = r_analytical(filename_r_RK4[i], filename_v_RK4[i], h[i])
        r_num = r_numerical(filename_r_RK4[i])

        Delta_max[i] = np.max(np.abs(r_exact - r_num))

    r_err = 0
    for i in range(1, len(h)):
        r_err += 0.25*np.log10(Delta_max[i]/Delta_max[i-1])/np.log10(h[i]/h[i-1])

    print("Error convergence rate with RK4: r_err = {}".format(r_err))


if error_convergence_rate_Euler:
    filename_r_Euler = ["Euler_r_1_1dt.txt", "Euler_r_1_01dt.txt", "Euler_r_1_001dt.txt", "Euler_r_1_0001dt.txt", "Euler_r_1_00001dt.txt"]
    filename_v_Euler = ["Euler_v_1_1dt.txt", "Euler_v_1_01dt.txt", "Euler_v_1_001dt.txt", "Euler_v_1_0001dt.txt", "Euler_v_1_00001dt.txt"]

    h = [1., 0.1, 0.01, 0.001, 0.0001]

    Delta_max = np.empty(len(h))
    for i in range(len(h)):
        r_exact = r_analytical(filename_r_Euler[i], filename_v_Euler[i], h[i])
        r_num = r_numerical(filename_r_Euler[i])

        Delta_max[i] = np.max(np.abs(r_exact - r_num))

    r_err = 0
    for i in range(1, len(h)):
        r_err += 0.25*np.log10(Delta_max[i]/Delta_max[i-1])/np.log10(h[i]/h[i-1])

    print("Error convergence rate with Euler: r_err = {}".format(r_err))


#filename_r = "RK4_r_1_0001dt.txt"
#filename_v = "RK4_v_1_0001dt.txt"
#h = 0.001

#r_analytical(filename_r, filename_v, h)
#r_numerical(filename_r)
