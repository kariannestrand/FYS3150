import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

n = np.linspace(10,100,10)
t = np.array([161, 729, 1715, 3079, 4861, 6995, 9530, 12498, 15746, 19624])

reg = linregress(np.log10(n),np.log10(t))
print(reg)

plt.title('Scaling of number of required transformations with matrix size N')
plt.plot(np.log10(n), np.log10(t), label = 'Slope = 2.07')
plt.legend()
plt.xlabel('log10(N)', size = 10)
plt.ylabel('log10(T)', size = 10)
plt.tight_layout()
plt.savefig('scaling.pdf')
#plt.show()

v_1 = np.array([0.0, 0.26287, 0.422533, 0.42533, 0.26287, 0, -0.26287, -0.42533, -0.42533, -0.26287, 0.0])
v_5 = np.array([0.0, -0.36180, -0.42533, -0.13820, 0.26287, 0.44721, 0.26287, -0.13820, -0.42533, -0.36180, 0.0])
v_7 = np.array([0.0, 0.13820, 0.26287, 0.36180, 0.42533, 0.44721, 0.42533, 0.36180, 0.26287, 0.13820, 0.0])

x_10 = np.linspace(0,1,11)

plt.title('Vector elements $v_i$ with corresponding positions $\hat{x}_i$')
plt.plot(x_10,v_1)
plt.plot(x_10,v_5)
plt.plot(x_10,v_7)
plt.xlabel('$\hat{x}_i$')
plt.ylabel('$v_i$')
plt.savefig('n10.pdf')
plt.show()

x_100 = np.linspace(0,1,100)
