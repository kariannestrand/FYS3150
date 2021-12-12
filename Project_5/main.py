import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa

h = 0.005
M = int(1.0/h + 1.0);
T = 0.008
dt = 2.5e-5
N = int(T/dt)

probability_deviation = True
probability_time_evolution = False

if probability_deviation:
    U_cube_0_T = pa.cx_cube()
    U_cube_0_T.load("0_v0.bin")
    U_cube_0 = np.transpose(U_cube_0_T)


    U_cube_10_T = pa.cx_cube()
    U_cube_10_T.load("10_v0.bin")
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
    plt.plot(t, abs(1 - p_0_array), ".", label = "$v_0 = 0$")
    plt.plot(t, abs(1 - p_10_array), ".", label = "$v_0 = 10^{10}$")
    #plt.title("Deviation from the normalized probability as a function of time", size = 22)
    plt.xlabel("Time", size = 22)
    plt.ylabel("|1 - $\sum_{i, j}\ p_{ij}$|", size = 22)
    plt.xticks(size = 20)
    plt.yticks(size = 20); plt.yscale("log")
    plt.legend(fontsize = 15)
    plt.show()




if probability_time_evolution:
    U_cube_T = pa.cx_cube()
    U_cube_T.load("T_0002.bin")
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
        plt.savefig('plot_test_3' + str(i)+ '.pdf')
        plt.show()
    for i in range(len(p)):
        plt.imshow(U_Re[i])
        plt.title(U_Re_title[i])
        plt.savefig('plot_test_2' + str(i)+ '.pdf')
        plt.show()
    for i in range(len(p)):
        plt.imshow(U_Im[i])
        plt.title(U_Im_title[i])
        plt.savefig('plot_test' + str(i)+ '.pdf')
        plt.show()
        

    '''
    x = 0.8
    y = np.linspace(0+h, 1-h, M-2)
    plt.plot(y, p_0002[int(x*(M-2)), :])
    plt.show()
    '''