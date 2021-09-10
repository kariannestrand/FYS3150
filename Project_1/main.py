import numpy as np
import matplotlib.pyplot as plt


n_list = [10, 100, 1000, 10000]

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
labels_10 = exact_10.keys()
x_exact_10 = exact_10[labels_10[0]]
u_10 = exact_10[labels_10[1]]

exact_100 = exact_list[1]
exact_100 = exact_100.keys()
x_exact_100 = exact_100[labels_100[0]]
u_100 = exact_100[labels_100[1]]

exact_1000 = exact_list[2]
exact_1000 = exact_100.keys()
x_exact_1000 = exact_10[labels_1000[0]]
u_1000 = exact_10[labels_1000[1]]

approx_10000 = approx_list[3]
labels_10000 = approx_10000.keys()
x_approx_10000 = approx_10000[labels_10000[0]]
u_10000 = approx_10000[labels_10000[1]]

approx_10 = approx_list[0]
labels_10 = approx_10.keys()
x_approx_10 = approx_10[labels_10[0]]
v_10 = approx_10[labels_10[1]]

approx_100 = approx_list[1]
labels_100 = approx_100.keys()
x_approx_100 = approx_10[labels_100[0]]
v_100 = approx_10[labels_100[1]]

approx_1000 = approx_list[2]
labels_1000 = approx_100.keys()
x_approx_1000 = approx_10[labels_100[0]]
v_1000 = approx_10[labels_100[1]]

approx_10000 = approx_list[3]
labels_10000 = approx_10000.keys()
x_approx_10000 = approx_10000[labels_10000[0]]
v_10000 = approx_10000[labels_10000[1]]


plt.figure()
plt.plot(x_exact_10, u_10, linestyle="-", marker="o", color='blue')
plt.title('Plot of $u(x)=1-(1-e^{-10})x-e^{-10x}$')
plt.xlabel('x-values')
plt.ylabel('u(x)')
plt.savefig("comparison.pdf")

plt.figure()
plt.plot(x_exact_10, u_10, linestyle="-", marker="o", color='blue')
plt.plot(x_approx_10, v_10, linestyle="-", marker="o", color='red')
plt.plot(x_approx_100, v_100, linestyle="-", color='green')
plt.plot(x_approx_1000, v_1000, linestyle="-", color='yellow')
plt.plot(x_approx_10000, v_10000, linestyle="-", color='brown')
plt.title('Comparison of exact and approximated solution')
plt.xlabel('x-values')
plt.ylabel('u(x)/v(x)')
plt.savefig("comparison.pdf")
#plt.show()