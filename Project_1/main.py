import numpy as np
import matplotlib.pyplot as plt

n_list = [10, 100, 1000, 10000, 100000, 1000000, 10000000]

def file_reader(filename):
    with open(filename, "r") as infile:
        labels = infile.readline().split()
        num_lists = len(labels)
        lines = infile.readlines()
        data = {key : [] for key in labels}
        for line in lines:
            vals = line.split()
            for i in range(num_lists):
                data[labels[i]].append(float(vals[i]))
    return data

exact_list = []
approx_list = []

for i in range(len(n_list)):
    filename = "exact_" + str(n_list[i]) + ".txt"

    data = file_reader(filename)
    labels = data.keys()
    data = {key : np.array(data[key]) for key in labels}
    exact_list.append(data)

for i in range(len(n_list)):
    filename = "approx_" + str(n_list[i]) + ".txt"

    data = file_reader(filename)
    labels = data.keys()
    data = {key : np.array(data[key]) for key in labels}
    approx_list.append(data)


exact_10 = exact_list[0]
labels_10 = list(exact_10.keys())
x_exact_10_meat = exact_10[labels_10[0]]
u_10_meat = exact_10[labels_10[1]]
x_exact_10_1 = np.append(x_exact_10_meat, 1)
x_exact_10 = np.insert(x_exact_10_1, 0, 0)
u_10_1 = np.append(u_10_meat, 0)
u_10 = np.insert(u_10_1, 0, 0)

exact_100 = exact_list[1]
labels_100 = list(exact_100.keys())
x_exact_100_meat = exact_100[labels_100[0]]
u_100_meat = exact_100[labels_100[1]]
x_exact_100_1 = np.append(x_exact_100_meat, 1)
x_exact_100 = np.insert(x_exact_100_1, 0, 0)
u_100_1 = np.append(u_100_meat, 0)
u_100 = np.insert(u_100_1, 0, 0)

exact_1000 = exact_list[2]
labels_1000 = list(exact_1000.keys())
x_exact_1000_meat = exact_1000[labels_1000[0]]
u_1000_meat = exact_1000[labels_1000[1]]
x_exact_1000_1 = np.append(x_exact_1000_meat, 1)
x_exact_1000 = np.insert(x_exact_1000_1, 0, 0)
u_1000_1 = np.append(u_1000_meat, 0)
u_1000 = np.insert(u_1000_1, 0, 0)

exact_10000 = exact_list[3]
labels_10000 = list(exact_10000.keys())
x_exact_10000_meat = exact_10000[labels_10000[0]]
u_10000_meat = exact_10000[labels_10000[1]]
x_exact_10000_1 = np.append(x_exact_10000_meat, 1)
x_exact_10000 = np.insert(x_exact_10000_1, 0, 0)
u_10000_1 = np.append(u_10000_meat, 0)
u_10000 = np.insert(u_10000_1, 0, 0)

exact_100000 = exact_list[4]
labels_100000 = list(exact_100000.keys())
x_exact_100000_meat = exact_100000[labels_100000[0]]
u_100000_meat = exact_100000[labels_100000[1]]
x_exact_100000_1 = np.append(x_exact_100000_meat, 1)
x_exact_100000 = np.insert(x_exact_100000_1, 0, 0)
u_100000_1 = np.append(u_100000_meat, 0)
u_100000 = np.insert(u_100000_1, 0, 0)

exact_1000000 = exact_list[5]
labels_1000000 = list(exact_1000000.keys())
x_exact_1000000_meat = exact_1000000[labels_1000000[0]]
u_1000000_meat = exact_1000000[labels_1000000[1]]
x_exact_1000000_1 = np.append(x_exact_1000000_meat, 1)
x_exact_1000000 = np.insert(x_exact_1000000_1, 0, 0)
u_1000000_1 = np.append(u_1000000_meat, 0)
u_1000000 = np.insert(u_1000000_1, 0, 0)

exact_10000000 = exact_list[6]
labels_10000000 = list(exact_10000000.keys())
x_exact_10000000_meat = exact_10000000[labels_10000000[0]]
u_10000000_meat = exact_10000000[labels_10000000[1]]
x_exact_10000000_1 = np.append(x_exact_10000000_meat, 1)
x_exact_10000000 = np.insert(x_exact_10000000_1, 0, 0)
u_10000000_1 = np.append(u_10000000_meat, 0)
u_10000000 = np.insert(u_10000000_1, 0, 0)

approx_10 = approx_list[0]
labels_10 = list(approx_10.keys())
x_approx_10_meat = approx_10[labels_10[0]]
v_10_meat = approx_10[labels_10[1]]
x_approx_10_1 = np.append(x_approx_10_meat, 1)
x_approx_10 = np.insert(x_approx_10_1, 0, 0)
v_10_1 = np.append(v_10_meat, 0)
v_10 = np.insert(v_10_1, 0, 0)

approx_100 = approx_list[1]
labels_100 = list(approx_100.keys())
x_approx_100_meat = approx_100[labels_100[0]]
v_100_meat = approx_100[labels_100[1]]
x_approx_100_1 = np.append(x_approx_100_meat, 1)
x_approx_100 = np.insert(x_approx_100_1, 0, 0)
v_100_1 = np.append(v_100_meat, 0)
v_100 = np.insert(v_100_1, 0, 0)


approx_1000 = approx_list[2]
labels_1000 = list(approx_1000.keys())
x_approx_1000_meat = approx_1000[labels_1000[0]]
v_1000_meat = approx_1000[labels_1000[1]]
x_approx_1000_1 = np.append(x_approx_1000_meat, 1)
x_approx_1000 = np.insert(x_approx_1000_1, 0, 0)
v_1000_1 = np.append(v_1000_meat, 0)
v_1000 = np.insert(v_1000_1, 0, 0)

approx_10000 = approx_list[3]
labels_10000 = list(approx_10000.keys())
x_approx_10000_meat = approx_10000[labels_10000[0]]
v_10000_meat = approx_10000[labels_10000[1]]
x_approx_10000_1 = np.append(x_approx_10000_meat, 1)
x_approx_10000 = np.insert(x_approx_10000_1, 0, 0)
v_10000_1 = np.append(v_10000_meat, 0)
v_10000 = np.insert(v_10000_1, 0, 0)

approx_100000 = approx_list[4]
labels_100000 = list(approx_100000.keys())
x_approx_100000_meat = approx_100000[labels_100000[0]]
v_100000_meat = approx_100000[labels_100000[1]]
x_approx_100000_1 = np.append(x_approx_100000_meat, 1)
x_approx_100000 = np.insert(x_approx_100000_1, 0, 0)
v_100000_1 = np.append(v_100000_meat, 0)
v_100000 = np.insert(v_100000_1, 0, 0)

approx_1000000 = approx_list[5]
labels_1000000 = list(approx_1000000.keys())
x_approx_1000000_meat = approx_1000000[labels_1000000[0]]
v_1000000_meat = approx_1000000[labels_1000000[1]]
x_approx_1000000_1 = np.append(x_approx_1000000_meat, 1)
x_approx_1000000 = np.insert(x_approx_1000000_1, 0, 0)
v_1000000_1 = np.append(v_1000000_meat, 0)
v_1000000 = np.insert(v_1000000_1, 0, 0)

approx_10000000 = approx_list[6]
labels_10000000 = list(approx_10000000.keys())
x_approx_10000000_meat = approx_10000000[labels_10000000[0]]
v_10000000_meat = approx_10000000[labels_10000000[1]]
x_approx_10000000_1 = np.append(x_approx_10000000_meat, 1)
x_approx_10000000 = np.insert(x_approx_10000000_1, 0, 0)
v_10000000_1 = np.append(v_10000000_meat, 0)
v_10000000 = np.insert(v_10000000_1, 0, 0)

abs_error_10 = np.abs(u_10-v_10)
abs_error_100 = np.abs(u_100-v_100)
abs_error_1000 = np.abs(u_1000-v_1000)
abs_error_10000 = np.abs(u_10000-v_10000)

rel_error_10 = np.abs((u_10_meat-v_10_meat)/u_10_meat)
rel_error_100 = np.abs((u_100_meat-v_100_meat)/u_100_meat)
rel_error_1000 = np.abs((u_1000_meat-v_1000_meat)/u_1000_meat)
rel_error_10000 = np.abs((u_10000_meat-v_10000_meat)/u_10000_meat)


rel_error_100000 = np.abs((u_100000_meat-v_100000_meat)/u_100000_meat)
rel_error_1000000 = np.abs((u_1000000_meat-v_1000000_meat)/u_1000000_meat)
rel_error_10000000 = np.abs((u_10000000_meat-v_10000000_meat)/u_10000000_meat)

log_abs_error_10 = np.log(abs_error_10)
log_abs_error_100 = np.log(abs_error_100)
log_abs_error_1000 = np.log(abs_error_1000)
log_abs_error_10000 = np.log(abs_error_10000)

log_rel_error_10 = np.log(rel_error_10)
log_rel_error_100 = np.log(rel_error_100)
log_rel_error_1000 = np.log(rel_error_1000)
log_rel_error_10000 = np.log(rel_error_10000)
log_rel_error_100000 = np.log(rel_error_100000)
log_rel_error_1000000 = np.log(rel_error_1000000)
log_rel_error_10000000 = np.log(rel_error_10000000)


max_rel_error_10 = np.amax(rel_error_10)
max_rel_error_100 = np.amax(rel_error_100)
max_rel_error_1000 = np.amax(rel_error_1000)
max_rel_error_10000 = np.amax(rel_error_10000)
max_rel_error_100000 = np.amax(rel_error_100000)
max_rel_error_1000000 = np.amax(rel_error_1000000)
max_rel_error_10000000 = np.amax(rel_error_10000000)

max_rel_error = np.array([max_rel_error_10, max_rel_error_100, max_rel_error_1000, max_rel_error_10000, max_rel_error_100000, max_rel_error_1000000, max_rel_error_10000000])
#print(max_rel_error)

n = np.linspace(10,10**7,7)

figwidth = 5.5
figheight = figwidth / 1.33333

plot = True
if plot:
    plt.figure(figsize=(figwidth, figheight))
    plt.plot(x_exact_10, u_10, linestyle="-", marker="o", color='blue')
    plt.title('Plot of exact solution $u(x)=1-(1-e^{-10})x-e^{-10x}$', fontsize=10)
    plt.xlabel('x-values')
    plt.ylabel('u(x)')
    plt.savefig("exact.pdf")
    plt.show()

    plt.figure(figsize=(figwidth, figheight))
    plt.plot(x_exact_10, u_10, '--', color='blue', label='Exact solution')
    plt.plot(x_approx_10, v_10, '--', color='black', label='Approximated solution n=10')
    plt.plot(x_approx_100, v_100, '--', color='green', label='Approximated solution n=10^2')
    plt.plot(x_approx_1000, v_1000, '--', color='orange', label='Approximated solution n=10^3')
    plt.plot(x_approx_10000, v_10000, '--', color='red', label='Approximated solution n=10^4')
    plt.title('Comparison of exact and approximated solutions to (write expression here)', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('u(x)/v(x)')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.savefig("comparison.pdf")
    plt.show()

    plt.figure(figsize=(figwidth, figheight))
    plt.plot(x_exact_10, abs_error_10, '--', color='blue', label='n=10')
    plt.plot(x_exact_100, abs_error_100, '--', color='red', label='n=10^2')
    plt.plot(x_exact_1000, abs_error_1000, '--', color='green', label='n=10^3')
    plt.plot(x_exact_10000, abs_error_10000, '--', color='yellow', label='n=10^4')
    plt.title("Absolute error for approximation of (write expression here)", fontsize=12)
    plt.xlabel("x")
    plt.ylabel("Absolute error")
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.savefig("abs_error.pdf")
    plt.show()

    plt.figure(figsize=(figwidth, figheight))
    plt.plot(x_exact_10_meat, rel_error_10, '--', color='blue', label='n=10')
    plt.plot(x_exact_100_meat, rel_error_100, '--', color='red', label='n=10^2')
    plt.plot(x_exact_1000_meat, rel_error_1000, '--', color='green', label='n=10^3')
    plt.plot(x_exact_10000_meat, rel_error_10000, '--', color='yellow', label='n=10^4')
    plt.title("Relative error", fontsize=12)
    plt.xlabel("x")
    plt.ylabel("Relative error for approximation of (write expression here)")
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.savefig("rel_error.pdf")
    plt.show()

    plt.figure(figsize=(figwidth, figheight))
    plt.plot(n, max_rel_error, '--')
    plt.xlabel('n')
    plt.ylabel('Max relative error')
    plt.title("Maximum relative error for approximation of (write expression here)", fontsize=12)
    plt.savefig("max_rel_error.pdf")
    plt.show()
