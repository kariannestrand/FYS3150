import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

burn_in_time = False
histogram = False
phase_transition = False
critical_temperature = False

# investigating the burnin time
if burn_in_time:
    filename1 = "eps_exp_24_ordered.bin"
    filename2 = "eps_exp_1_ordered.bin"
    filename3 = "eps_exp_24_unordered.bin"
    filename4 = "eps_exp_1_unordered.bin"
    filename5 = "m_exp_1_ordered.bin"
    filename6 = "m_exp_1_unordered.bin"
    filename7 = "m_exp_24_ordered.bin"
    filename8 = "m_exp_24_unordered.bin"

    eps_exp24_ord = np.loadtxt(filename1, skiprows = 0)
    eps_exp1_ord = np.loadtxt(filename2, skiprows = 0)
    eps_exp24_un = np.loadtxt(filename3, skiprows = 0)
    eps_exp1_un = np.loadtxt(filename4, skiprows = 0)
    m_exp1_ord = np.loadtxt(filename5, skiprows = 0)
    m_exp1_un = np.loadtxt(filename6, skiprows = 0)
    m_exp24_ord = np.loadtxt(filename7, skiprows = 0)
    m_exp24_un = np.loadtxt(filename8, skiprows = 0)

    n = len(eps_exp24_ord)  #number of cycles
    T = np.linspace(0, n, n)

    plt.plot(T, eps_exp24_ord, label = "T = 2.4 J/$k_B$ - Ordered")
    plt.plot(T, eps_exp24_un, label = "T = 2.4 J/$k_B$ - Unordered")
    plt.plot(T, eps_exp1_ord, label = "T = 1.0 J/$k_B$ - Ordered")
    plt.plot(T, eps_exp1_un, label = "T = 1.0 J/$k_B$ - Unordered")
    plt.ylabel("Expected energy per spin/[J]", size = 20)
    plt.xlabel("Number of Monte Carlo cycles", size = 20)
    plt.xticks(size = 17); plt.xscale("log")
    plt.yticks(size = 17)
    plt.legend(fontsize = 14)
    plt.show()


    plt.plot(T, m_exp24_ord, label = "T = 2.4 J/$k_B$ - Ordered")
    plt.plot(T, m_exp24_un, label = "T = 2.4 J/$k_B$ - Unordered")
    plt.plot(T, m_exp1_ord, label = "T = 1.0 J/$k_B$ - Ordered")
    plt.plot(T, m_exp1_un, label = "T = 1.0 J/$k_B$ - Unordered")
    plt.ylabel("Expected absolute magnetization per spin", size = 17)
    plt.xlabel("Number of Monte Carlo cycles", size = 20)
    plt.xticks(size = 17); plt.xscale("log")
    plt.yticks(size = 17)
    plt.legend(fontsize = 14)
    plt.show()

# making histograms of samples of the expectation value of the energy
if histogram:
    filename1 = "eps_histo_1_unord.bin"
    filename24 = "eps_histo_24_unord.bin"

    eps_1 = np.loadtxt(filename1, skiprows = 0)
    eps_24 = np.loadtxt(filename24, skiprows = 0)

    n, bins, patches = plt.hist(eps_1, 120, density = True)
    plt.title("Probability of energy per spin with T = 1.0 J/$k_B$", size = 22)
    plt.ylabel('Normalized frequency', fontsize = 20)
    plt.xlabel('Energy per spin/[J]', fontsize = 20)
    plt.xticks(fontsize = 17)
    plt.yticks(fontsize = 17)
    plt.legend()
    plt.show()

    n, bins, patches = plt.hist(eps_24, 120, density = True)
    plt.title("Probability of energy per spin with T = 2.4 J/$k_B$", size = 22)
    plt.ylabel('Normalized frequency', fontsize = 20)
    plt.xlabel('Energy per spin/[J]', fontsize = 20)
    plt.xticks(fontsize = 17)
    plt.yticks(fontsize = 17)
    plt.legend()
    plt.show()

# investigating phace transitions
if phase_transition:
    L = ["40", "60", "80", "100"]

    filename_list_eps = ["eps_exp_40L.bin", "eps_exp_60L.bin", "eps_exp_80L.bin", "eps_exp_100L.bin"]
    for i in range(len(filename_list_eps)):
        filename = filename_list_eps[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        eps_exp_unsorted = loadtxt[:, 0]
        T_unsorted = loadtxt[:, 1]
        T, e_exp = zip(*sorted(zip(T_unsorted, eps_exp_unsorted)))
        plt.plot(T, e_exp, '.', label = "L = " + L[i])
        plt.xlabel("Temperature/[J/$k_B$]", size = 20)
        plt.ylabel("Expected energy per spin/[J]", size = 20)
        plt.xticks(size = 17)
        plt.yticks(size = 17)
        plt.legend(fontsize = 14)
    plt.show()


    filename_list_m = ["m_exp_40L_noburnin.bin", "m_exp_60L_noburnin.bin", "m_exp_80L_noburnin.bin", "m_exp_100L_noburnin.bin"]
    for i in range(len(filename_list_m)):
        filename = filename_list_m[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        m_exp_unsorted = loadtxt[:, 0]
        T_unsorted = loadtxt[:, 1]
        T, m_exp = zip(*sorted(zip(T_unsorted, m_exp_unsorted)))
        plt.plot(T, m_exp, '.', label = "L = " + L[i])
        plt.xlabel("Temperature/[J/$k_B$]", size = 20)
        plt.ylabel("Expected absolute magnetization per spin", size = 20)
        plt.xticks(size = 17)
        plt.yticks(size = 17)
        plt.legend(fontsize = 14)
    plt.show()


    filename_list_cV = ["cV_40L.bin", "cV_60L.bin", "cV_80L.bin", "cV_100L.bin"]
    for i in range(len(filename_list_cV)):
        filename = filename_list_cV[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        cV_unsorted = loadtxt[:, 0]
        T_unsorted = loadtxt[:, 1]
        T, cV = zip(*sorted(zip(T_unsorted, cV_unsorted)))
        plt.plot(T, cV, '.', label = "L = " + L[i])
        plt.xlabel("Temperature/[J/$k_B$]", size = 20)
        plt.ylabel("Heat capacity per spin/[$k_B$]", size = 20)
        plt.xticks(size = 17)
        plt.yticks(size = 17)
        plt.legend(fontsize = 14)
    plt.show()



    filename_list_chi = ["chi_40L_noburnin.bin", "chi_60L_noburnin.bin", "chi_80L_noburnin.bin", "chi_100L_noburnin.bin"]
    for i in range(len(filename_list_chi)):
        filename = filename_list_chi[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        chi_unsorted = loadtxt[:, 0]
        T_unsorted = loadtxt[:, 1]
        T, chi = zip(*sorted(zip(T_unsorted, chi_unsorted)))
        plt.plot(T, chi, '.', label = "L = " + L[i])
        plt.xlabel("Temperature/[J/$k_B$]", size = 20)
        plt.ylabel("Suceptibility per spin/[1/J]", size = 20)
        plt.xticks(size = 17)
        plt.yticks(size = 17)
        plt.legend(fontsize = 14)
    plt.show()

# estimating critical temperature
if critical_temperature:
    T_max_list = []

    L = ["40", "60", "80", "100"]
    filename_list_chi = ["chi_40L.bin", "chi_60L.bin", "chi_80L.bin", "chi_100L.bin"]
    for i in range(len(filename_list_chi)):
        filename = filename_list_chi[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        chi_unsorted = loadtxt[:, 0]
        T_unsorted = loadtxt[:, 1]
        T, chi = zip(*sorted(zip(T_unsorted, chi_unsorted)))
        chi_max = np.amax(chi)
        chi_max_index = np.argmax(chi)
        T_max = T[chi_max_index]
        T_max_list.append(T_max)
        plt.plot(T, chi)
        plt.plot(T_max, chi_max, "o", label = "$T_c$(" + L[i] + ") = " + str(T_max) + " J/$k_B$")
        plt.xlabel("Temperature/[J/$k_B$]", size = 20)
        plt.ylabel("Suceptibility per spin/[1/J]", size = 20)
        plt.xticks(size = 17)
        plt.yticks(size = 17)
        plt.legend(fontsize = 17)
    plt.show()

    L_array = np.array([40, 60, 80, 100])
    T_max_array = np.array(T_max_list)
    for i in range(len(T_max_list)):
        plt.plot(1./float(L[i]), T_max_list[i], "o", label = "")

    x = 1./L_array
    y = T_max_array
    m, b = np.polyfit(x, y, 1)
    slope, intercept, r, p, se = linregress(x, y)
    Tc = "{0:.3f}".format(intercept)
    plt.plot(x, m*x + b, label = "T$_c$ = " + str(Tc))
    plt.title("Estimate of T$_c$ using suceptibility", size = 22)
    plt.xlabel("Inverce lattice size", size = 20)
    plt.ylabel("Critical temperature/[J/$k_B$]", size = 20)
    plt.xticks(size = 17)
    plt.yticks(size = 17)
    plt.legend(fontsize = 17)
    plt.show()
