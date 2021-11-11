import numpy as np
import matplotlib.pyplot as plt

phace_transition = True

# investigating phace transitions
if phace_transition:
    filename = "eps_exp.bin"
    loadtxt = np.loadtxt(filename, skiprows = 1)
    eps_exp = loadtxt[:, 0]
    T = loadtxt[:, 1]
    plt.plot(T, eps_exp)
    plt.show()


    filename = "m_exp.bin"
    loadtxt = np.loadtxt(filename, skiprows = 1)
    m_exp = loadtxt[:, 0]
    T = loadtxt[:, 1]
    plt.plot(T, m_exp)
    plt.show()


    filename = "cV.bin"
    loadtxt = np.loadtxt(filename, skiprows = 1)
    cV = loadtxt[:, 0]
    T = loadtxt[:, 1]
    plt.plot(T, cV)
    plt.show()

    filename = "chi.bin"
    loadtxt = np.loadtxt(filename, skiprows = 1)
    chi = loadtxt[:, 0]
    T = loadtxt[:, 1]
    plt.plot(T, chi)
    plt.show()
