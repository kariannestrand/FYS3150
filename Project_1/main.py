import numpy as np
import matplotlib.pyplot as plt


n = [10, 10**2, 10**3, 10**4, 10**5]


""" plotting exact solution for all n if True """
exact_solution = False
if exact_solution:

    # reading from file, extracting x and u from exact_100.txt:
    filename = "exact_1000.txt"
    exact = np.loadtxt(filename, skiprows = 1)
    x = exact[:, 0]
    u = exact[:, 1]

    # plotting exact solution using x and u:
    plt.plot(x, u, label = "n = 1000")
    plt.title("Exact solution", fontsize = 20)
    plt.xlabel("x", size = 20); plt.xticks(size = 15)
    plt.ylabel("u(x)", size = 20); plt.yticks(size = 15)
    plt.legend(fontsize = 15)
    plt.tight_layout()
    plt.savefig("exact_solution.pdf")
    plt.show()



""" plotting comparison between exact and approximate solutions for all n if True """
compare_solution = True
if compare_solution:

    # reading from files, extracting u from exact_n.txt:
    filename_exact = "exact_1000.txt"
    exact = np.loadtxt(filename_exact, skiprows = 1)
    u = exact[:, 1]

    for i in range(len(n)):
        # reading from files, extracting x and v from approx_n.txt:
        filename_approx = "approx_" + str(n[i]) + ".txt"
        approx = np.loadtxt(filename_approx, skiprows = 1)
        x = approx[:, 0]
        v = approx[:, 1]

        # plotting compared solutions using x, u and v:
        if n[i] == 1000:
            plt.plot(x, u, label = "u(x), n = 1000")
        plt.plot(x, v, "--", label = "v(x), n = " + str(n[i]))
        plt.title("Comparison between exact and approximated solutions", fontsize = 14)
        plt.xlabel("x", size = 20); plt.xticks(size = 15)
        plt.ylabel("u(x) and v(x)", size = 20); plt.yticks(size = 15)
        plt.legend(fontsize = 10)
        plt.tight_layout()
        plt.savefig("compare_solution.pdf")
    plt.show()



""" plotting absolute error for all n if True """
absolute_error = False
if absolute_error:

    # function calculating absolute error:
    def abs_error(u, v):
        return np.log(np.abs(u - v))

    for i in range(len(n)):
        # reading from files, extracting x, u from exact_n.txt and v from approx_n.txt:
        filename_exact = "exact_" + str(n[i]) + ".txt"
        exact = np.loadtxt(filename_exact, skiprows = 1)
        x = exact[:, 0]
        u = exact[:, 1]

        filename_approx = "approx_" + str(n[i]) + ".txt"
        approx = np.loadtxt(filename_approx, skiprows = 1)
        v = approx[:, 1]

        # plotting absolute error by calling function abs_error with u and v:
        plt.plot(x, abs_error(u, v), "--", label = "n = " + str(n[i]))
        plt.title("Absolute error", fontsize = 20)
        plt.xlabel("x", size = 20); plt.xticks(size = 15)
        plt.ylabel("$log_{10}(\Delta)$", size = 20); plt.yticks(size = 15)
        plt.tight_layout()
        plt.legend(fontsize = 15)
        plt.savefig("absolute_error.pdf")
    plt.show()



""" plotting relative error for all n if True """
relative_error = False
if relative_error:

    # function calculating relative error:
    def rel_error(u, v):
        return np.log(np.abs((u - v)/u))

    for i in range(len(n)):
        # reading from files, extracting x, u from exact_n.txt and v from approx_n.txt:
        filename_exact = "exact_" + str(n[i]) + ".txt"
        exact = np.loadtxt(filename_exact, skiprows = 1)
        x = exact[:, 0]
        u = exact[:, 1]

        filename_approx = "approx_" + str(n[i]) + ".txt"
        approx = np.loadtxt(filename_approx, skiprows = 1)
        v = approx[:, 1]

        # plotting relative error by calling function rel_error with u and v:
        plt.plot(x, rel_error(u, v), "--", label = "n = " + str(n[i]))
        plt.title("Relative error", fontsize = 20)
        plt.xlabel("x", size = 20); plt.xticks(size = 15)
        plt.ylabel("$log_{10}(\epsilon)$", size = 20); plt.yticks(size = 15)
        plt.legend(loc = "upper left", fontsize = 10)
        plt.tight_layout()
        plt.savefig("relative_error.pdf")
    plt.show()



""" plotting maxiumum relative error for n up to 10**7 if True """
maximum_relative_error = False
if maximum_relative_error:

    # function calculating maximum relative error:
    def max_rel_error(u, v):
        return np.amax(np.log(np.abs((u - v)/u)))

    # defining n up to 10**7:
    n = [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7]


    for i in range(len(n)):
        # reading from files, extracting u from exact_n.txt and v from approx_n.txt:
        filename_exact = "exact_" + str(n[i]) + ".txt"
        exact = np.loadtxt(filename_exact, skiprows = 1)
        u = exact[:, 1]

        filename_approx = "approx_" + str(n[i]) + ".txt"
        approx = np.loadtxt(filename_approx, skiprows = 1)
        v = approx[:, 1]

        # plotting maximum relative error as a function of n by calling function max_rel_error with u and v:
        plt.plot(n[i], max_rel_error(u, v), "o", label = "n = " + str(n[i]))
        plt.title("Maximum relative error", fontsize = 20)
        plt.xlabel("n", size = 20); plt.xticks(size = 15)
        plt.ylabel("$log_{10}(\epsilon)_{max}$", size = 20); plt.yticks(size = 15)
        plt.xscale("log")
        plt.tight_layout()
        plt.legend(fontsize = 15)
        plt.savefig("maximum_relative_error.pdf")
    plt.show()
