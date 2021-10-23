import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pyarma as pa

# constants
q = 1.
m = 40.078

B0 = 9.65e1
V0 = 9.65e8
d = 1e4

omega_0 = q*B0/m
omega_z = np.sqrt(2*q*V0/(m*d**2))

omega_p = 0.5*(omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2))
omega_m = 0.5*(omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2))

# loading files for dt = 0.002
V_1 = pa.mat()
V_1.load("Euler_v_0_0002dt.bin")
vel_1 = np.transpose(V_1)

vy_1 = vel_1[:, 1]

R_1 = pa.mat()
R_1.load("Euler_r_0_0002dt.bin")
pos_1 = np.transpose(R_1)

x_1 = pos_1[:, 0]
y_1 = pos_1[:, 1]
z_1 = pos_1[:, 2]

t_1 = np.linspace(0,100,len(x_1))

#analytic solution dt = 0.002

N_1 = len(x_1)

x_0_1 = x_1[0]
v_0_1 = vy_1[0]
z_0_1 = z_1[0]

Ap_1 = (v_0_1 + omega_m*x_0_1)/(omega_m - omega_p)
Am_1 = - (v_0_1 + omega_p*x_0_1)/(omega_m - omega_p)

x_ex_1 = np.empty(N_1)
y_ex_1 = np.empty(N_1)
z_ex_1 = np.empty(N_1)

r_exact_1 = np.empty(N_1)
h_1 = 0.002

for i in range(N_1):
    x_ex_1[i] = Ap_1*np.cos(omega_p*i*h_1) + Am_1*np.cos(omega_m*i*h_1)
    y_ex_1[i] = - Ap_1*np.sin(omega_p*i*h_1) - Am_1*np.sin(omega_m*i*h_1)
    z_ex_1[i] = z_0_1*np.cos(omega_z*i*h_1)


r_1 = np.sqrt((x_ex_1-x_1)**2+(y_ex_1-y_1)**2+(z_ex_1-z_1)**2)
r_1_exact = np.sqrt(x_ex_1[i]**2 + y_ex_1[i]**2 + z_ex_1[i]**2)

rel_err_1 = r_1/r_1_exact


# loading files for dt = 0.001
V_2 = pa.mat()
V_2.load("Euler_v_0_0001dt.bin")
vel_2 = np.transpose(V_2)

vy_2 = vel_2[:, 1]

R_2 = pa.mat()
R_2.load("Euler_r_0_0001dt.bin")
pos_2 = np.transpose(R_2)

x_2 = pos_2[:, 0]
y_2 = pos_2[:, 1]
z_2 = pos_2[:, 2]

t_2 = np.linspace(0,100,len(x_2))

#analytic solution dt = 0.001

N_2 = len(x_2)

x_0_2 = x_2[0]
v_0_2 = vy_2[0]
z_0_2 = z_2[0]

Ap_2 = (v_0_2 + omega_m*x_0_2)/(omega_m - omega_p)
Am_2 = - (v_0_2 + omega_p*x_0_2)/(omega_m - omega_p)

x_ex_2 = np.empty(N_2)
y_ex_2 = np.empty(N_2)
z_ex_2 = np.empty(N_2)

h_2 = 0.001

for i in range(N_2):
    x_ex_2[i] = Ap_2*np.cos(omega_p*i*h_2) + Am_2*np.cos(omega_m*i*h_2)
    y_ex_2[i] = - Ap_2*np.sin(omega_p*i*h_2) - Am_2*np.sin(omega_m*i*h_2)
    z_ex_2[i] = z_0_2*np.cos(omega_z*i*h_2)

r_2 = np.sqrt((x_ex_2-x_2)**2+(y_ex_2-y_2)**2+(z_ex_2-z_2)**2)
r_2_exact = np.sqrt(x_ex_2[i]**2 + y_ex_2[i]**2 + z_ex_2[i]**2)

rel_err_2 = r_2/r_2_exact


# loading files for dt = 0.0008
V_3 = pa.mat()
V_3.load("Euler_v_0_00008dt.bin")
vel_3 = np.transpose(V_3)

vy_3 = vel_3[:, 1]

R_3 = pa.mat()
R_3.load("Euler_r_0_00008dt.bin")
pos_3 = np.transpose(R_3)

x_3 = pos_3[:, 0]
y_3 = pos_3[:, 1]
z_3 = pos_3[:, 2]

t_3 = np.linspace(0,100,len(x_3))

#analytic solution dt = 0.0008

N_3 = len(x_3)

x_0_3 = x_3[0]
v_0_3 = vy_3[0]
z_0_3 = z_3[0]

Ap_3 = (v_0_3 + omega_m*x_0_3)/(omega_m - omega_p)
Am_3 = - (v_0_3 + omega_p*x_0_3)/(omega_m - omega_p)

x_ex_3 = np.empty(N_3)
y_ex_3 = np.empty(N_3)
z_ex_3 = np.empty(N_3)

h_3 = 0.0008

for i in range(N_3):
    x_ex_3[i] = Ap_3*np.cos(omega_p*i*h_3) + Am_3*np.cos(omega_m*i*h_3)
    y_ex_3[i] = - Ap_3*np.sin(omega_p*i*h_3) - Am_3*np.sin(omega_m*i*h_3)
    z_ex_3[i] = z_0_3*np.cos(omega_z*i*h_3)

r_3 = np.sqrt((x_ex_3-x_3)**2+(y_ex_3-y_3)**2+(z_ex_3-z_3)**2)
r_3_exact = np.sqrt(x_ex_3[i]**2 + y_ex_3[i]**2 + z_ex_3[i]**2)

rel_err_3 = r_3/r_3_exact


# loading files for dt = 0.0004
V_4 = pa.mat()
V_4.load("Euler_v_0_00004dt.bin")
vel_4 = np.transpose(V_4)

vy_4 = vel_4[:, 1]

R_4 = pa.mat()
R_4.load("Euler_r_0_00004dt.bin")
pos_4 = np.transpose(R_4)

x_4 = pos_4[:, 0]
y_4 = pos_4[:, 1]
z_4 = pos_4[:, 2]

t_4 = np.linspace(0,100,len(x_4))

#analytic solution dt = 0.0004

N_4 = len(x_4)

x_0_4 = x_4[0]
v_0_4 = vy_4[0]
z_0_4 = z_4[0]

Ap_4 = (v_0_4 + omega_m*x_0_4)/(omega_m - omega_p)
Am_4 = - (v_0_4 + omega_p*x_0_4)/(omega_m - omega_p)

x_ex_4 = np.empty(N_4)
y_ex_4 = np.empty(N_4)
z_ex_4 = np.empty(N_4)

h_4 = 0.0004

for i in range(N_4):
    x_ex_4[i] = Ap_4*np.cos(omega_p*i*h_4) + Am_4*np.cos(omega_m*i*h_4)
    y_ex_4[i] = - Ap_4*np.sin(omega_p*i*h_4) - Am_4*np.sin(omega_m*i*h_4)
    z_ex_4[i] = z_0_4*np.cos(omega_z*i*h_4)


r_4 = np.sqrt((x_ex_4-x_4)**2+(y_ex_4-y_4)**2+(z_ex_4-z_4)**2)
r_4_exact = np.sqrt(x_ex_4[i]**2 + y_ex_4[i]**2 + z_ex_4[i]**2)

rel_err_4 = r_4/r_4_exact



# loading files for dt = 0.0001
V_5 = pa.mat()
V_5.load("Euler_v_0_00001dt.bin")
vel_5 = np.transpose(V_5)

vy_5 = vel_5[:, 1]

R_5 = pa.mat()
R_5.load("Euler_r_0_00001dt.bin")
pos_5 = np.transpose(R_5)

x_5 = pos_5[:, 0]
y_5 = pos_5[:, 1]
z_5 = pos_5[:, 2]

t_5 = np.linspace(0,100,len(x_5))

#analytic solution dt = 0.0001

N_5 = len(x_5)

x_0_5 = x_5[0]
v_0_5 = vy_5[0]
z_0_5 = z_5[0]

Ap_5 = (v_0_5 + omega_m*x_0_5)/(omega_m - omega_p)
Am_5 = - (v_0_5 + omega_p*x_0_5)/(omega_m - omega_p)

x_ex_5 = np.empty(N_5)
y_ex_5 = np.empty(N_5)
z_ex_5 = np.empty(N_5)

h_5 = 0.0001

for i in range(N_5):
    x_ex_5[i] = Ap_5*np.cos(omega_p*i*h_5) + Am_5*np.cos(omega_m*i*h_5)
    y_ex_5[i] = - Ap_5*np.sin(omega_p*i*h_5) - Am_5*np.sin(omega_m*i*h_5)
    z_ex_5[i] = z_0_5*np.cos(omega_z*i*h_5)

r_5 = np.sqrt((x_ex_5-x_5)**2+(y_ex_5-y_5)**2+(z_ex_5-z_5)**2)
r_5_exact = np.sqrt(x_ex_5[i]**2 + y_ex_5[i]**2 + z_ex_5[i]**2)

rel_err_5 = r_5/r_5_exact

plt.title('Relative error with Euler')
plt.plot(t_1, rel_err_1, label = 'dt = 0.002')
plt.plot(t_2, rel_err_2, label = 'dt = 0.001')
plt.plot(t_3, rel_err_3, label = 'dt = 0.0008')
plt.plot(t_4, rel_err_4, label = 'dt = 0.0004')
plt.plot(t_5, rel_err_5, label = 'dt = 0.0001')
plt.xlabel('t[$\mu s$]')
plt.ylabel('$\epsilon$')
plt.legend()
plt.savefig('rel_error_euler.pdf')

delta_max_1 = np.max(r_1)
delta_max_2 = np.max(r_2)
delta_max_3 = np.max(r_3)
delta_max_4 = np.max(r_4)
delta_max_5 = np.max(r_5)

delta_max = [delta_max_1, delta_max_2, delta_max_3, delta_max_4, delta_max_5]
dt_list = [0.015, 0.008, 0.004, 0.001, 0.0008]

sum = 0
for k in range(1,5):
    sum += (np.log10(delta_max[k]/delta_max[k-1]))/(np.log10(float(dt_list[k])/float(dt_list[k-1])))
r_err = (1/4)*sum

print(r_err)


