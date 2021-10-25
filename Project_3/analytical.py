import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pyarma as pa

# constants
q = 1.
m = 40.078

B0 = 9.65e1
V0 = 9.65e8
d = 1e4

omega_0 = q*B0/m
omega_z = np.sqrt(2*q*V0/(m*d**2))

omega_p = 0.5*(omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2))
omega_m = 0.5*(omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2))

analytical = True

relative_error_RK4 = False
relative_error_Euler = False

error_convergence_rate_RK4 = True
error_convergence_rate_Euler = True

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
    x_ex = np.empty(N)
    y_ex = np.empty(N)
    z_ex = np.empty(N)

    r_exact = np.empty(N)

    for i in range(N):
        x_ex[i] = Ap*np.cos(omega_p*i*h) + Am*np.cos(omega_m*i*h)
        y_ex[i] = - Ap*np.sin(omega_p*i*h) - Am*np.sin(omega_m*i*h)
        z_ex[i] = z_0*np.cos(omega_z*i*h)

        r_exact[i] = np.sqrt(x_ex[i]**2 + y_ex[i]**2 + z_ex[i]**2)

    return r_exact, x_ex, y_ex, z_ex

def relative_error(filename_r, filename_v, h):

    r = np.sqrt((r_analytical(filename_r, filename_v, h)[1]-R(filename_r)[0])**2 + (r_analytical(filename_r, filename_v, h)[2]-R(filename_r)[1])**2 + (r_analytical(filename_r, filename_v, h)[3]-R(filename_r)[2])**2)
    r_exact = r_analytical(filename_r, filename_v, h)[0]

    rel_err = r/r_exact

    return rel_err, r

def error_convergence_rate(filename_r, filename_v, h):
    N = len(h)
    delta_max = np.empty(N)

    for i in range(N):
        r = relative_error(filename_r[i], filename_v[i], h[i])[1]
        delta_max[i] = np.max(r)

    print(delta_max)

    sum = 0
    for k in range(1,5):
        sum += (np.log10(delta_max[k]/delta_max[k-1]))/(np.log10(float(dt_list[k])/float(dt_list[k-1])))
    r_err = (1/4)*sum

    return r_err

def time(filename_r):
    N = len(R(filename_r)[0])
    t = np.linspace(0, 100, N)
    return t



if analytical:
    filename_r = 'z1.bin'
    filename_v = 'z2.bin'
    h = 0.001
    t = time(filename_r)
    Rz = R(filename_r)[2]
    z_ex = r_analytical(filename_r, filename_v, h)[3]
    plt.plot(t, z_ex, label = 'Analytical solution')
    plt.plot(t, Rz, '--', label = 'Numerical solution')
    plt.xlabel('z/[$\mu$m]')
    plt.ylabel('t/[$\mu$s]')
    plt.legend()
    plt.savefig('zt_1.pdf')
    plt.show()

if relative_error_RK4:
    filename_r_list = ["RK4_r_0_0016dt.bin", "RK4_r_0_0008dt.bin", "RK4_r_0_0004dt.bin", "RK4_r_0_0002dt.bin", "RK4_r_0_0001dt.bin"]
    filename_v_list = ["RK4_v_0_0016dt.bin", "RK4_v_0_0008dt.bin", "RK4_v_0_0004dt.bin", "RK4_v_0_0002dt.bin", "RK4_v_0_0001dt.bin"]
    dt_list = [0.016, 0.008, 0.004, 0.002, 0.001]

    N = len(dt_list)
    for i in range(N):
        t = time(filename_r_list[i])
        rel_err = relative_error(filename_r_list[i], filename_v_list[i], dt_list[i])[0]
        plt.plot(t, rel_err, label = "h = " + str(dt_list[i]))

    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("$\epsilon$", size = 12)
    plt.legend()
    plt.savefig('rel_error_rk4.pdf')

if relative_error_Euler:
    filename_r_list = ["Euler_r_0_00016dt.bin", "Euler_r_0_00008dt.bin", "Euler_r_0_00004dt.bin", "Euler_r_0_00002dt.bin", "Euler_r_0_00001dt.bin"]
    filename_v_list = ["Euler_v_0_00016dt.bin", "Euler_v_0_00008dt.bin", "Euler_v_0_00004dt.bin", "Euler_v_0_00002dt.bin", "Euler_v_0_00001dt.bin"]
    dt_list = [0.0016, 0.0008, 0.0004, 0.0002, 0.0001]


    N = len(dt_list)
    for i in range(N):
        t = time(filename_r_list[i])
        rel_err = relative_error(filename_r_list[i], filename_v_list[i], dt_list[i])[0]
        plt.plot(t, rel_err, label = "h = " + str(dt_list[i]))

    plt.xlabel("t/[$\mu$s]", size = 12)
    plt.ylabel("$\epsilon$", size = 12)
    plt.legend()
    plt.savefig('rel_error_euler.pdf')


if error_convergence_rate_RK4:
    filename_r_list = ["RK4_r_0_0016dt.bin", "RK4_r_0_0008dt.bin", "RK4_r_0_0004dt.bin", "RK4_r_0_0002dt.bin", "RK4_r_0_0001dt.bin"]
    filename_v_list = ["RK4_v_0_0016dt.bin", "RK4_v_0_0008dt.bin", "RK4_v_0_0004dt.bin", "RK4_v_0_0002dt.bin", "RK4_v_0_0001dt.bin"]
    dt_list = [0.016, 0.008, 0.004, 0.002, 0.001]

    r_err = error_convergence_rate(filename_r_list, filename_v_list, dt_list)
    print("Error convergence rate with RK4: r_err = {}".format(r_err))

if error_convergence_rate_Euler:
    filename_r_list = ["Euler_r_0_00016dt.bin", "Euler_r_0_00008dt.bin", "Euler_r_0_00004dt.bin", "Euler_r_0_00002dt.bin", "Euler_r_0_00001dt.bin"]
    filename_v_list = ["Euler_v_0_00016dt.bin", "Euler_v_0_00008dt.bin", "Euler_v_0_00004dt.bin", "Euler_v_0_00002dt.bin", "Euler_v_0_00001dt.bin"]
    dt_list = [0.0016, 0.0008, 0.0004, 0.0002, 0.0001]

    r_err = error_convergence_rate(filename_r_list, filename_v_list, dt_list)
    print("Error convergence rate with Euler: r_err = {}".format(r_err))
