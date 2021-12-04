import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

N = np.linspace(10,100,10) # Size of matrix
t = np.array([161, 729, 1715, 3079, 4861, 6995, 9530, 12498, 15746, 19624]) # Number of required transformations

reg = linregress(np.log10(N),np.log10(t)) # Linear regression to find the scaling exponent

# Plotting the scaling of number of required transformations with matrix size N

plt.title('Scaling of number of required transformations with matrix size N')
plt.plot(np.log10(N), np.log10(t), label = 'Slope = 2.07')
plt.legend()
plt.xlabel('log10(N)', size = 10)
plt.ylabel('log10(T)', size = 10)
plt.tight_layout()
plt.savefig('scaling.pdf')
#plt.show()

# Plotting the vector elements v with corresponding positions x for n = 100

x_10 = np.linspace(0,1,11) # Positions x_i

# Vector elements from vectors corresponding to the lowest eigenvalues

v_1_10 = np.array([0.0, 0.26287, 0.422533, 0.42533, 0.26287, 0, -0.26287, -0.42533, -0.42533, -0.26287, 0.0])
v_2_10 = np.array([0.0, -0.36180, -0.42533, -0.13820, 0.26287, 0.44721, 0.26287, -0.13820, -0.42533, -0.36180, 0.0])
v_3_10 = np.array([0.0, 0.13820, 0.26287, 0.36180, 0.42533, 0.44721, 0.42533, 0.36180, 0.26287, 0.13820, 0.0])

plt.title('Vector elements $v_i$ with corresponding positions $\hat{x}_i$')
plt.plot(x_10,v_1_10)
plt.plot(x_10,v_2_10)
plt.plot(x_10,v_3_10)
plt.xlabel('$\hat{x}_i$')
plt.ylabel('$v_i$')
plt.savefig('n10.pdf')
#plt.show()

# Plotting the vector elements v with corresponding positions x for n = 100

x_100 = np.linspace(0,1,101) # Positions x_i

# Loading vector elements from file

v_1 = np.loadtxt('eigvec_1.txt')
v_2 = np.loadtxt('eigvec_2.txt')
v_3 = np.loadtxt('eigvec_3.txt')

# Adding the solution v_0=0 and v_n=0 to arrays

v_11 = np.insert(v_1, 0, 0)
v_21 = np.insert(v_2, 0, 0)
v_31 = np.insert(v_3, 0, 0)
v_1_100 = np.insert(v_11, 100, 0)
v_2_100 = np.insert(v_21, 100, 0)
v_3_100 = np.insert(v_31, 100, 0)

plt.title('Vector elements $v_i$ with corresponding positions $\hat{x}_i$')
plt.plot(x_100,v_1_100)
plt.plot(x_100,v_2_100)
plt.plot(x_100,v_3_100)
plt.xlabel('$\hat{x}_i$')
plt.ylabel('$v_i$')
plt.savefig('n100.pdf')
#plt.show()
