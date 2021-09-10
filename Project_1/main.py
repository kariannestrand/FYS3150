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

abs_error_10 = np.abs(u_10-v_10)
abs_error_100 = np.abs(u_100-v_100)
abs_error_1000 = np.abs(u_1000-v_1000)
abs_error_10000 = np.abs(u_10000-v_10000)

rel_error_10 = np.abs((u_10_meat-v_10_meat)/u_10_meat)
rel_error_100 = np.abs((u_100_meat-v_100_meat)/u_100_meat)
rel_error_1000 = np.abs((u_1000_meat-v_1000_meat)/u_1000_meat)
rel_error_10000 = np.abs((u_10000_meat-v_10000_meat)/u_10000_meat)

log_abs_error_10 = np.log(abs_error_10)
log_abs_error_100 = np.log(abs_error_100)
log_abs_error_1000 = np.log(abs_error_1000)
log_abs_error_10000 = np.log(abs_error_10000)

log_rel_error_10 = np.log(rel_error_10)
log_rel_error_100 = np.log(rel_error_100)
log_rel_error_1000 = np.log(rel_error_1000)
log_rel_error_10000 = np.log(rel_error_10000)

max_rel_error_10 = np.amax(rel_error_10)
max_rel_error_100 = np.amax(rel_error_100)
max_rel_error_1000 = np.amax(rel_error_1000)
max_rel_error_10000 = np.amax(rel_error_10000)

max_rel_error = np.array([max_rel_error_10, max_rel_error_100, max_rel_error_1000, max_rel_error_10000])

n = np.linspace(1,4,4)
print(n)


plt.plot(x_exact_10, u_10, linestyle="-", marker="o", color='blue')
plt.title('Plot of $u(x)=1-(1-e^{-10})x-e^{-10x}$')
plt.xlabel('x-values')
plt.ylabel('u(x)')
plt.savefig("exact.pdf")
plt.show()

plt.plot(x_exact_10, u_10, linestyle="-", marker="o", color='blue')
plt.plot(x_approx_10, v_10, linestyle="-", marker="o", color='red')
plt.plot(x_approx_100, v_100, linestyle="-", color='green')
plt.plot(x_approx_1000, v_1000, linestyle="-", color='yellow')
plt.plot(x_approx_10000, v_10000, linestyle="-", color='brown')
plt.title('Comparison of exact and approximated solution')
plt.xlabel('x')
plt.ylabel('u(x)/v(x)')
plt.savefig("comparison.pdf")
plt.show()

plt.plot(x_exact_10, log_abs_error_10, color='blue')
plt.plot(x_exact_100, log_abs_error_100, color='red')
plt.plot(x_exact_1000, log_abs_error_1000, color='green')
plt.plot(x_exact_10000, log_abs_error_10000, color='yellow')
plt.savefig("abs_error.pdf")
plt.show()

plt.plot(x_exact_10_meat, log_rel_error_10, color='blue')
plt.plot(x_exact_100_meat, log_rel_error_100, color='red')
plt.plot(x_exact_1000_meat, log_rel_error_1000, color='green')
plt.plot(x_exact_10000_meat, log_rel_error_10000, color='yellow')
plt.savefig("rel_error.pdf")
plt.show()

plt.plot(n, max_rel_error)
plt.show()
