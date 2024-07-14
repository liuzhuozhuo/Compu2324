import numpy as np
from numba import jit
import pandas as pd
import matplotlib.pyplot as plt
import timeit

machine = "mac"
start = timeit.default_timer() 

N = 25
L = np.sqrt(N)
h = 0.002
epsilon = 1
sigma = 1
m = 1
v_0 = 0
k=1

skip = 50 # Cada cuantos pasos se guarda la informacion

Time = 600
time = np.arange(0, Time, h)

def init_cond():
    r = np.zeros((N, 2))
    v =np.zeros((N, 2))
    for i in range(N):
        r[i] = np.array([i%L+0.5, i//L +0.5])
        theta = np.random.rand()*2*np.pi
        v[i] = v_0 * np.array([np.sin(theta), np.cos(theta)])
    return r, v

def cond_contorno(r):
    if(r[0] > L):
        r[0] = r[0] - L
    if(r[0] < 0):
        r[0] = r[0] + L
    if(r[1] > L):
        r[1] = r[1] - L
    if(r[1] < 0):
        r[1] = r[1] + L
    return r

def cond_contorno_distancia(r):
    if(r[0] > L/2):
        r[0] = r[0] - L
    elif(r[0] < -L/2):
        r[0] = r[0] + L
    if(r[1] > L/2):
        r[1] = r[1] - L
    elif(r[1] < -L/2):
        r[1] = r[1] + L
    return r

def lennard_jones(r):
    R = compute_distance(r)
    acc = np.zeros((N, 2))
    for i in range(N):
        for j in range(N):
            if(i!=j):
                norm  = np.linalg.norm(R[i, j])
                if (norm < 3):
                    acc[i] = acc[i] + 4*R[i, j]* epsilon * (6*np.power((sigma/norm), 5) - 12*np.power((sigma/norm), 11))/(norm*m)
    return acc, R

def compute_distance(r):
    R = np.zeros((N, N, 2))
    for i in range(0, N-1):
        for j in range(i+1, N):
            R[i, j] = r[j]- r[i]
            R[i, j] = cond_contorno_distancia(R[i, j])
            R[j, i] = -R[i, j]
    return R

def verlet_algorithm(r, v, a):
    w = np.zeros((N, 2))
    for i in range(N):
        r[i] = r[i] + h*v[i] + h*h*a[i]/2
        r[i] = cond_contorno(r[i])
        w[i] = v[i] + h*a[i]/2
    a, R = lennard_jones(r)
    for i in range(N):
        v[i] = w[i] + h*a[i]/2
    return r, v, a, R

def compute_energy(v, R):
    T = 0
    V = 0
    for i in range(N):
        T = T + 0.5*m*np.linalg.norm(v[i])**2
        for j in range(N):
            if(i!=j):
                norm  = np.linalg.norm(R[i, j])
                V = V + 4*epsilon * (np.power((sigma/norm), 12) - np.power((sigma/norm), 6))
    return T, V

def compute_average_speed(v):
    v_prom = 0
    for i in range(N):
        v_prom = v_prom + np.linalg.norm(v[i])
    return v_prom/N

def write_vector(r, f):
    for i in range(N):
        f.write(f"{r[i, 0]}, {r[i, 1]}\n")

def write_velocity(v, f):
    for i in range(N):
        f.write(f"{v[i, 0]}, {v[i, 1]}, {np.linalg.norm(v[i])}\n")

def compute_mean_square_displacement(r1, r2):
    num = 0
    if len(r1) != len(r2):
        for i in range(len(r1)):
            num = np.linalg.norm(r1[i] - r2)**2
    else:
        for i in range(len(r1)):
            num = np.linalg.norm(r1[i] - r2[i])**2
    return num/(len(r1))

f = open(f"voluntario1/optimizaciones/no_optimizado/posiciones_{v_0}.txt", "w")
f_energia = open(f"voluntario1/optimizaciones/no_optimizado/energias_{v_0}.txt", "w")

T, V = np.zeros(len(time)), np.zeros(len(time))

r, v = init_cond()
write_vector(r, f)
a, R = lennard_jones(r)
T[0], V[0] = compute_energy(v, R)
f_energia.write(f"{T[0]}, {V[0]}\n")

for i in range(1, len(time)):
    r, v, a, R = verlet_algorithm(r, v, a)
    T[i], V[i] = compute_energy(v, R)
    if (time[i]%60 == 0):
        v = v*1.1
    if (i%skip == 0):
        f.write(f"\n")
        write_vector(r, f)
        f_energia.write(f"{T[i]}, {V[i]}\n")

f.close()
f_energia.close()

file_in = f"voluntario1/optimizaciones/no_optimizado/posiciones_{v_0}.txt"

# Lectura del fichero de datos
# ========================================
# Lee el fichero a una cadena de texto
with open(file_in, "r") as f:
    data_str = f.read()

# Inicializa la lista con los datos de cada fotograma.
# frames_data[j] contiene los datos del fotograma j-ésimo
frames_data = list()

# Itera sobre los bloques de texto separados por líneas vacías
# (cada bloque corresponde a un instante de tiempo)
for frame_data_str in data_str.split("\n\n"):
    # Inicializa la lista con la posición de cada planeta
    frame_data = list()

    # Itera sobre las líneas del bloque
    # (cada línea da la posición de un planta)
    for planet_pos_str in frame_data_str.split("\n"):
        # Lee la componente x e y de la línea
        planet_pos = np.fromstring(planet_pos_str, sep=",")
        # Si la línea no está vacía, añade planet_pos a la lista de 
        # posiciones del fotograma
        if planet_pos.size > 0:
            frame_data.append(np.fromstring(planet_pos_str, sep=","))

    # Añade los datos de este fotograma a la lista
    frames_data.append(frame_data)

data = np.array(frames_data)

f_desplazamiento = open(f"voluntario1/optimizaciones/no_optimizado/desplazamiento_{v_0}.txt", "w")

time_between_measurement = 0.01
desplazamiento = np.zeros(int(Time/(time_between_measurement*skip)))
n = int(len(data[:,0])/len(desplazamiento))
for i in range (len(desplazamiento)):
    desplazamiento[i] = compute_mean_square_displacement(data[:,0][i*n:(i+1)*n], data[:,5][i*n:(i+1)*n])
    f_desplazamiento.write(f"{desplazamiento[i]}\n")

f_desplazamiento.close()

stop = timeit.default_timer()
timer = stop - start

f_timer = open(f"voluntario1/optimizaciones/no_optimizado/datos_{machine}.txt", "a+")
f_timer.write(f"{N}, {h}, {Time}, {timer}\n")
f_timer.close()

