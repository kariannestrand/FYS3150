import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# one particle, no interactions

T_1 = pa.mat()
T_1.load("r_0_0001.bin")
R_1 = np.transpose(T_1)

t = np.linspace(0, 100, len(R_1))

x_1 = R_1[:, 0]
y_1 = R_1[:, 1]
z_1 = R_1[:, 2]

#print(R_1)

#plt.plot(x_1, y_1)
#plt.plot(t, z_1)
#plt.savefig('test.pdf')

ax = plt.axes(projection="3d")
plt.tight_layout()
ax.plot3D(x_1, y_1, z_1, 'blue', label='Trajectory of particle 1')
#plt.savefig('3d_test.pdf')
plt.show()



'''
# two interacting particles
T_1 = pa.mat()
T_2 = pa.mat()
T_1.load("r_0_0001.bin")
T_2.load("r_1_0001.bin")
R_1 = np.transpose(T_1)
R_2 = np.transpose(T_2)

S_1 = pa.mat()
S_2 = pa.mat()
S_1.load("v_0_0001.bin")
S_2.load("v_1_0001.bin")
V_1 = np.transpose(S_1)
V_2 = np.transpose(S_2)


t = np.linspace(0,100,len(R_1))

x_1 = R_1[:, 0]
y_1 = R_1[:, 1]
z_1 = R_1[:, 2]

x_2 = R_2[:, 0]
y_2 = R_2[:, 1]
z_2 = R_2[:, 2]


#plt.plot(x_1, y_1)
#plt.plot(x_2, y_2)
plt.plot(t, z_1)
plt.savefig('test.pdf')


'''
