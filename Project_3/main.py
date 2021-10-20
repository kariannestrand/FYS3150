import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

n = 1
interaction = False
read_velocity = False
read_position = False

z_t = False
x_y = False

v_x = False
v_y = False
v_z = False

trajectory = False

analytical = False

relative_error_RK4 = True
relative_error_Euler = True



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
    t = np.linspace(0,100,len(Rz1))
    fig = plt.figure()
    plt.plot(t, Rz1)

    if n == 1:
        plt.title("Movement of one particle in the z-direction", size = 10)

    if n == 2:
        plt.plot(t, Rz2)

        if interaction:
            plt.title("Movement of two particles in the z-direction w/ interaction", size = 10)
        else:
            plt.title("Movement of two particles in the z-direction w/o interaction", size = 10)

    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("z/[$\mu$m]", size = 12)

    plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

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


    if n == 1:
        plt.title("Trajectory of one particle", size = 12)

    if n == 2:
        ax.plot3D(Rx2, Ry2, Rz2, 'red', label ='Trajectory of particle 2')

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



def r_analytical(filename_r, filename_v):

    pos = np.loadtxt(filename_r, skiprows = 1)
    Rx1 = pos[:, 0]
    Rz1 = pos[:, 2]

    vel = np.loadtxt(filename_v, skiprows = 1)
    Vy1 = pos[:, 1]

    x_0 = Rx1[0]
    v_0 = Vy1[0]
    z_0 = Rz1[0]

    q = 1.
    m = 40.078

    B0 = 9.65e1
    V0 = 9.65e8
    d = 10e4

    omega_0 = q*B0/m
    omega_z = np.sqrt(2*q*V0/(m*d**2))

    omega_p = 0.5*(omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2));
    omega_m = 0.5*(omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2));

    Ap = (v_0 + omega_p*x_0)/(omega_p - omega_m)
    Am = - (v_0 + omega_m*x_0)/(omega_p - omega_m)


    # analytical solution, r_exact
    r_exact_list = []
    for i in range(len(Rx1)):
        dt = 0.001
        x = Ap*np.cos(omega_p*i*dt) + Am*np.cos(omega_m*i*dt)
        y = - Ap*np.sin(omega_p*i*dt) - Am*np.sin(omega_m*i*dt)

        z_Re = z_0 + np.cos(omega_z*i*dt)
        z_Im = np.sin(omega_z*dt)

        # append length of r, or vector r (x, y, z) ?? anf z_Re or what??
        r_exact_list.append(np.sqrt(x**2 + y**2 + z_Re**2))

    r_exact = np.array(r_exact_list)

    return r_exact

def r_numerical(filename_r):
    pos = np.loadtxt(filename_r, skiprows = 1)
    Rx1 = pos[:, 0]
    Ry1 = pos[:, 1]
    Rz1 = pos[:, 2]

    r_num = np.empty(len(Rx1))
    for i in range(len(Rx1)):
        r_num[i] = np.sqrt(Rx1[i]**2 + Ry1[i]**2 + Rz1[i]**2)

    return r_num

def relative_error(filename_r, filename_v):
    r_exact = r_analytical(filename_r, filename_v)
    r_num = r_numerical(filename_r)

    rel_err = np.abs((r_exact - r_num)/r_num)

    return rel_err

def time(filename_r):
    N = len(r_numerical(filename_r))
    t = np.linspace(0, 100, N)
    return t


if relative_error_RK4:
    filename_r_RK4 = ["RK4_r_1_1dt.txt", "RK4_r_1_01dt.txt", "RK4_r_1_001dt.txt", "RK4_r_1_0001dt.txt", "RK4_r_1_00001dt.txt"]
    filename_v_RK4 = ["RK4_v_1_1dt.txt", "RK4_v_1_01dt.txt", "RK4_v_1_001dt.txt", "RK4_v_1_0001dt.txt", "RK4_v_1_00001dt.txt"]


    for i in range(len(filename_r_RK4)):
        plt.plot(time(filename_r_RK4[i]), relative_error(filename_r_RK4[i], filename_v_RK4[i]))
    plt.show()


if relative_error_Euler:
    filename_r_Euler = ["Euler_r_1_1dt.txt", "Euler_r_1_01dt.txt", "Euler_r_1_001dt.txt", "Euler_r_1_0001dt.txt", "Euler_r_1_00001dt.txt"]
    filename_v_Euler = ["Euler_v_1_1dt.txt", "Euler_v_1_01dt.txt", "Euler_v_1_001dt.txt", "Euler_v_1_0001dt.txt", "Euler_v_1_00001dt.txt"]

    for i in range(len(filename_r_Euler)):
        plt.plot(time(filename_r_Euler[i]), relative_error(filename_r_Euler[i], filename_v_Euler[i]))
    plt.show()




    """
    Delta_max_list = []
    dt = [1., 0.1, 0.01, 0.001, 0.0001]
    for j in range(len(dt)):

        Delta_list = []
        relative_error_list = []

        N = len(Rx1)
        for i in range(N):
            x = Ap*np.cos(omega_p*i*dt[j]) + Am*np.cos(omega_m*i*dt[j])
            y = - Ap*np.sin(omega_p*i*dt[j]) - Am*np.sin(omega_m*i*dt[j])

            z_Re = z_0 + np.cos(omega_z*i*dt[j])
            z_Im = np.sin(omega_z*dt[j])

            # append length of r, or vector r (x, y, z) ??
            r_exact = np.sqrt(x**2 + y**2 + z_Re**2)
            r_num = np.sqrt(Rx1[i]**2 + Ry1[i]**2 + Rz1[i]**2)

            relative_error_list.append(np.abs((r_exact - r_num)/r_num))
            Delta_list.append(np.abs(r_exact - r_num))

        Delta_max_list.append(np.amax(Delta_list))


    Delta_max = np.array(Delta_max_list)

    r_err = 0
    for i in range(1, len(dt)):
        r_err += 0.25*np.log10(Delta_max[i]/Delta_max[i-1])/np.log10(dt[i]/dt[i-1])
    print(r_err)
    """




"""
for i in range(n_particles):
    filename = "RK4_r.txt"
    Rx = []
    Ry = []
    Rz = []
    with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split()
            if len(particles) > 1:
                Rx.append(particles[i])
"""
