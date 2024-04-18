import numpy as np
import matplotlib.pyplot as plt

#Inicializa las posiciones iniciales aleatorios
def init_spins (N):
    spins = np.random.choice([-1, 1], (N,N))
    return spins

#Función que cambia el valor del spin
def moves_spin (spins, N, T):
    #generate a random position
    pos = np.random.randint(0, N, size=(1, 2))[0]
    d_energy = 2*spins[pos[0], pos[1]]*(spins[(pos[0]+1)%N, pos[1]] + spins[(pos[0]-1)%N, pos[1]] + spins[(pos[0]), (pos[1]-1)%N] + spins[(pos[0]), (pos[1]+1)%N] )
    p = min(1, np.exp(-d_energy/T))
    xi = np.random.rand()
    if(xi < p):
        spins[pos[0], pos[1]] = -spins[pos[0], pos[1]]
    return spins

#Función que devuelve la energía del sistema
def energy(spins, N):
    sum = 0
    for i in range(N):
        for j in range(N):
            sum += spins[i, j]*(spins[(i+1)%N, j] + spins[(i-1)%N, j] + spins[(i), (j-1)%N] + spins[(i), (j+1)%N])
    return -sum/2

#Función que devuelve la probabilidad 
def probability (energy, T, k):
    return np.exp((-energy)/(k*T))

#Función que devuelve la magnetización del sistema
def magnetization (spins, N):
    sum = 0
    for i in range(N):
        for j in range(N):
            sum += spins[i, j]
    return sum/(N**2)

#Función que devuelve el promedio de la magnetización
def total_magnetization (T, N, time, dt, k):
    step = dt * N**2
    n_data = int(time/dt)

    spins = init_spins(N)
    all_spins = np.zeros((n_data, N, N))
    prob = np.zeros(n_data)
    m = np.zeros(n_data)

    z = 0
    M = 0

    i =0
    while (i < n_data):
        j = 0
        mag = 0
        while (j < step):
            spins = moves_spin(spins, N, T)
            mag += magnetization(spins, N)
            j+=1
        all_spins[i] = spins
        ener = energy(spins, N)
        m[i] = mag/step
        prob[i]=probability(ener, T, k)
        M += m[i]*prob[i]
        z += prob[i]
        i+=1

    return M/z

N = 8
T = np.arange(0.5, 5, 0.5)
k=1
time = 10**6
dt = 100

mag = np.zeros(len(T))

i = 0
while(i < len(T)):
    mag[i] = total_magnetization(T[i], N, time, dt, k)
    i+=1

fig=plt.figure(figsize=(7,7)) #Size of the plot
ax=fig.add_subplot(111)
for i in range(len(T)):
    plt.scatter(T[i], mag[i], label = f'{T[i]}')
plt.legend()
plt.grid()

fig.savefig("{}.png".format('obligatoria2/magnetizacion'))