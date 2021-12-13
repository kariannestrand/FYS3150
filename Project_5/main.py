from matplotlib.colors import Normalize
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa

h = 0.005
M = int(1.0/h + 1.0);
T = 0.002
dt = 2.5e-5
N = int(T/dt)

probability_deviation =  False
probability_time_evolution = True

if probability_deviation:
    U_cube_0_T = pa.cx_cube()
    U_cube_0_T.load("no_barrier.bin")
    U_cube_0 = np.transpose(U_cube_0_T)


    U_cube_10_T = pa.cx_cube()
    U_cube_10_T.load("double_slit.bin")
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
    #plt.title("Deviation from the normalized probability as a function of time", size = 22)
    plt.xlabel("Time", size = 24)
    plt.ylabel("|1 - $\sum_{i, j}\ p_{ij}$|", size = 22)
    plt.xticks([0.0, 0.002, 0.004, 0.006, 0.008], size = 24)
    plt.yticks(size = 24); plt.yscale("log")
    plt.legend(fontsize = 24)
    plt.tight_layout()
    plt.savefig('probability.pdf')
    plt.show()




if probability_time_evolution:
    U_cube_T = pa.cx_cube()
    U_cube_T.load("colourmap.bin")
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


    #saved x and y from c++ to file to make the grid from 0 to 1
    x = pa.mat()
    y = pa.mat()
    x.load("x.bin")
    y.load("y.bin")
    X, Y = np.meshgrid(x,y)
    
    for i in range(len(p)):
        plt.contourf(X, Y, p[i], 20)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xticks([0.2,0.4,0.6,0.8,1])
        plt.yticks([0,0.2,0.4,0.6,0.8])
        cb = plt.colorbar()
        cb.set_label(label='Probability')
        cb.ax.tick_params()
        plt.savefig('prob' + str(i)+ '.pdf')
        plt.show()
        plt.close("all")
    
    
    for i in range(len(p)):
        plt.contourf(X, Y, U_Re[i], 20)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xticks([0.2,0.4,0.6,0.8,1])
        plt.yticks([0,0.2,0.4,0.6,0.8])
        cb = plt.colorbar()
        cb.set_label(label='Real part')
        cb.ax.tick_params()
        plt.savefig('real' + str(i)+ '.pdf')
        plt.show()
        plt.close("all")

    for i in range(len(p)):
        plt.contourf(X, Y, U_Im[i], 20)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xticks([0.2,0.4,0.6,0.8,1])
        plt.yticks([0,0.2,0.4,0.6,0.8])
        cb = plt.colorbar()
        cb.set_label(label='Imaginary part')
        cb.ax.tick_params()
        plt.savefig('imag' + str(i)+ '.pdf')
        plt.show()
        plt.close("all")



'''
this is what we had before, saved it here just in case
    for i in range(len(p)):
        plt.ax.imshow(U_Im[i])
        plt.title(U_Im_title[i])
        plt.colorbar()
        plt.savefig('plot_test' + str(i)+ '.pdf')
        plt.show()
'''
    
    
    