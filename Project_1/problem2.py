import numpy as np
import matplotlib.pyplot as plt


n_list = [10, 100, 1000, 10000]

def file_reader(n_list):
    data = np.zeros((len(4), 2)) 
    for i in range(len(n_list)):
        v = np.zeros(n_list[i]+1)
        filename = "problem7_" + str(n_list[i]) + ".txt"
        v = 

x_exact = []; u_exact = []
with open("problem2.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x_exact.append(float(line[0])); u_exact.append(float(line[1]))

x10 = []; u10 = []
with open("problem7_10.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x10.append(float(line[0])); u10.append(float(line[1]))

x100 = []; u100 = []
with open("problem7_100.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x100.append(float(line[0])); u100.append(float(line[1]))

x1000 = []; u1000 = []
with open("problem7_1000.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x1000.append(float(line[0])); u1000.append(float(line[1]))

x10000 = []; u10000 = []
with open("problem7_10000.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x10000.append(float(line[0])); u10000.append(float(line[1]))


plt.plot(x_exact, u_exact, linestyle="-", marker="o", color='blue')
plt.plot(x10, u10, linestyle="-", marker="o", color='red')
plt.plot(x100, u100, linestyle="-", color='green')
plt.plot(x1000, u1000, linestyle="-", color='yellow')
plt.plot(x10000, u10000, linestyle="-", color='brown')
plt.title('Plot of $u(x)=1-(1-e^{-10})x-e^{-10x}$')
plt.xlabel('x-values')
plt.ylabel('u(x)')
plt.savefig("problem7.pdf")
plt.show()