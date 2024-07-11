import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import timeit

machine = "windows"

start = timeit.default_timer() 

g = 9.81
l = 1
h = 0.001

time = 100000

# Usado para las graficas
dot = '\u0307'

E = np.array([1, 3, 5, 10, 15])

phi_init = np.array([0, np.pi/16, np.pi/8, np.pi/6, np.pi/4, np.pi/3])
psi_init = np.array([0, np.pi/16, np.pi/8, np.pi/4, np.pi/3, np.pi/2])

n_interation = int(time/h)

@jit(nopython = True)
def polar_to_cart (phi):
    x = np.zeros((2, len(phi)))
    for i in range(len(phi)):
        x[0, i] = l * np.cos(phi[i]-np.pi/2)
        x[1, i] = l * np.sin(phi[i]-np.pi/2)
    return x

@jit(nopython = True)
def momento_to_velocidad (vec):
    div = 2-(np.cos(vec[0]-vec[1]))**2
    phi_dot = (vec[2]-vec[3]*np.cos(vec[0]-vec[1]))/div
    psi_dot = (2*vec[3]-vec[2]*np.cos(vec[0]-vec[1]))/div
    return np.array([vec[0], vec[1], phi_dot, psi_dot])

@jit(nopython = True)
def function (vec):
    div = 2-(np.cos(vec[0]-vec[1]))**2
    pphi_dot = 2*np.sin(vec[1]-vec[0])*(vec[2]*vec[3]*(np.cos(vec[0]-vec[1]))**2 - (2*vec[3]**2 + vec[2]**2)*np.cos(vec[0]-vec[1]) + 2*vec[2]*vec[3])/div**2 - 2*g*np.sin(vec[0])
    ppsi_dot = 2*np.sin(vec[0]-vec[1])*(vec[2]*vec[3]*(np.cos(vec[0]-vec[1]))**2 - (2*vec[3]**2 + vec[2]**2)*np.cos(vec[0]-vec[1]) + 2*vec[2]*vec[3])/div**2 - g*np.sin(vec[1]) 
    phi_dot = (vec[2]-vec[3]*np.cos(vec[0]-vec[1]))/div
    psi_dot = (2*vec[3]-vec[2]*np.cos(vec[0]-vec[1]))/div
    return np.array([phi_dot, psi_dot, pphi_dot, ppsi_dot])

@jit(nopython = True)
def Runge_Kutta (vec_angulos, h):
    k_1 = h*function(vec_angulos)
    k_2 = h*function(vec_angulos + k_1/2)
    k_3 = h*function(vec_angulos + k_2/2)
    k_4 = h*function(vec_angulos + k_3)
    return vec_angulos + (k_1 + k_2 + k_3 + k_4)/6

@jit(nopython = True)
def run_code(n, j, k):
    vec_Cohete_total = np.zeros((int(n_interation), 4))

    phi_dot_init = np.sqrt(E[n] - 2*g*(1-np.cos(phi_init[j]))- g*(1-np.cos(psi_init[k])))
    ppsi_init = phi_dot_init*(np.cos(psi_init[k]-phi_init[j]))**2
    pphi_init = 2*phi_dot_init

    vec_Cohete_total[0] = np.array([phi_init[j], psi_init[k], pphi_init, ppsi_init])
    for i in range(n_interation-1):
        vec_Cohete_total[i+1] = Runge_Kutta(vec_Cohete_total[i], h)
    return vec_Cohete_total 

vec_pendulo = run_code(4, 0, 0)
pos_pendulo_1 = polar_to_cart(vec_pendulo[:, 0])
pos_pendulo_2 = pos_pendulo_1 + polar_to_cart(vec_pendulo[:, 1])

f = open("voluntario2\datos_optimizacion\datos.txt", "w")
for i in range(int(n_interation)):
    if i%10 == 0:
        f.write(f"{pos_pendulo_1[0, i]}, {pos_pendulo_1[1, i]}\n")
        f.write(f"{pos_pendulo_2[0, i]}, {pos_pendulo_2[1, i]}\n")
        f.write(f"\n")
f.close()

stop = timeit.default_timer()
timer = stop - start

f_timer = open(f"voluntario2\datos_optimizacion\datos_{machine}.txt", "a+")
f_timer.write(f"{h}, {n_interation}, {time},  {timer}\n")
