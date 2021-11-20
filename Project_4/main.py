import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

burn_in_time = False
histogram = False
phase_transition = False

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

    loadtxt1 = np.loadtxt(filename1, skiprows = 0)
    loadtxt2 = np.loadtxt(filename2, skiprows = 0)
    loadtxt3 = np.loadtxt(filename3, skiprows = 0)
    loadtxt4 = np.loadtxt(filename4, skiprows = 0)
    loadtxt5 = np.loadtxt(filename5, skiprows = 0)
    loadtxt6 = np.loadtxt(filename6, skiprows = 0)
    loadtxt7 = np.loadtxt(filename7, skiprows = 0)
    loadtxt8 = np.loadtxt(filename8, skiprows = 0)
    
    eps_exp24_ord = loadtxt1
    eps_exp1_ord = loadtxt2
    eps_exp24_un = loadtxt3
    eps_exp1_un = loadtxt4
    m_exp1_ord = loadtxt5
    m_exp1_un = loadtxt6
    m_exp24_ord = loadtxt7
    m_exp24_un = loadtxt8
    n = 100000  #number of cycles
    T = np.linspace(0,1000001,n)/1000
    
    plt.plot(T, eps_exp24_ord)
    plt.plot(T, eps_exp24_un)
    plt.rc('xtick', labelsize=10) 
    plt.rc('ytick', labelsize=10)
    plt.ylabel('$<\epsilon>$') 
    plt.xlabel('t/[$10^3$ cycles]')
    plt.savefig('eps_burnin.pdf')
    plt.show()
    
    plt.plot(T, eps_exp1_ord)
    plt.plot(T, eps_exp1_un)
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    plt.show()
    plt.plot(T, m_exp1_ord)
    plt.plot(T, m_exp1_un)
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    plt.show()
    plt.plot(T, m_exp24_ord)
    plt.plot(T, m_exp24_un)
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    plt.show()
    
    


# making histograms of samples of the expectation value of the energy
if histogram:
    filename1 = 'eps_histo_1.bin'
    filename2 = 'eps_histo_24.bin'
    loadtxt1 = np.loadtxt(filename1, skiprows = 0)
    loadtxt2 = np.loadtxt(filename2, skiprows = 0)
    eps_1 = loadtxt1
    eps_2 = loadtxt2

    n, bins, patches = plt.hist(eps_1, 180, density=True)
    plt.xlim(-2.01,-1.8)
    plt.ylabel('Probability', fontsize=15)
    plt.xlabel('$\epsilon$', fontsize=15)
    plt.xticks(fontsize=10) 
    plt.yticks(fontsize=10)
    plt.savefig('histo_1.pdf')
    plt.show()
    n, bins, patches = plt.hist(eps_2, 120, density=True)
    plt.ylabel('Probability', fontsize=10)
    plt.xlabel('$\epsilon$', fontsize=10)
    plt.xticks(fontsize=10) 
    plt.yticks(fontsize=10)
    plt.savefig('histo_2.pdf')
    plt.show()

# for testing the files
filename = 'cV_60L_noburnin.bin'
loadtxt = np.loadtxt(filename, skiprows = 0)
eps = loadtxt[:, 0]
T = loadtxt[:, 1]
T, eps = zip(*sorted(zip(T, eps)))
T_max = T[np.argmax(eps)]
print(T_max)
plt.plot(T, eps, 'o')
plt.savefig('test.pdf')
plt.show()


# investigating phace transitions
if phase_transition:
    filename_list_eps = ["eps_exp_40L.bin", "eps_exp_60L.bin", "eps_exp_80L.bin", "eps_exp_100L.bin"]
    for i in range(len(filename_list_eps)):
        filename = filename_list_eps[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        eps_exp = loadtxt[:, 0]
        T = loadtxt[:, 1]
        T, cV = zip(*sorted(zip(T, eps_exp)))
        plt.plot(T, eps_exp, 'o')
        plt.show()


    filename_list_m = ["m_exp_40L.bin", "m_exp_60L.bin", "m_exp_80L.bin", "m_exp_100L.bin"]
    for i in range(len(filename_list_m)):
        filename = filename_list_m[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        m_exp = loadtxt[:, 0]
        T = loadtxt[:, 1]
        T, m_exp = zip(*sorted(zip(T, m_exp)))
        plt.plot(T, m_exp, 'o')
        plt.show()

    
    filename_list_cv = ["cV_40L.bin", "cV_60L.bin", "cV_80L.bin", "cV_100L.bin"]
    for i in range(len(filename_list_cv)):
        filename = filename_list_cv[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        cV = loadtxt[:, 0]
        T = loadtxt[:, 1]
        T, cV = zip(*sorted(zip(T, cV)))
        T_max = T[np.argmax(cV)]
        print(T_max)
        plt.plot(T, cV, 'o')
        plt.show()

    

    filename_list_chi = ["chi_40L.bin", "chi_60L.bin", "chi_80L.bin", "chi_100L.bin"]
    for i in range(len(filename_list_chi)):
        filename = filename_list_chi[i]
        loadtxt = np.loadtxt(filename, skiprows = 0)
        chi = loadtxt[:, 0]
        T = loadtxt[:, 1]
        T, chi = zip(*sorted(zip(T, chi)))
        T_max = T[np.argmax(chi)]
        print(T_max)
        plt.plot(T, chi, 'o')
        plt.show()
