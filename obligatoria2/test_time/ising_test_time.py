
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
import timeit

start = timeit.default_timer() 

N = 8
T = np.arange(0.5, 5, 0.5)
k = 1
time_ising = 10**4
dt_ising = 1
n_data = int(time_ising/dt_ising)
ising_data = np.zeros((len(T), n_data, N, N))

@jit(nopython=True)
def init_spins (N):
    spins = np.random.choice(np.array([-1, 1]), (N,N))
    return spins

#Función que cambia el valor del spin
@jit(nopython=True)
def moves_spin (spins, N, T):
    #generate a random position
    pos = np.random.randint(0, N, size=(1, 2))[0]
    d_energy = 2*spins[pos[0], pos[1]]*(spins[(pos[0]+1)%N, pos[1]] + spins[(pos[0]-1)%N, pos[1]] + spins[(pos[0]), (pos[1]-1)%N] + spins[(pos[0]), (pos[1]+1)%N] )
    p = min(1, np.exp(-d_energy/T))
    xi = np.random.rand()
    if(xi < p):
        spins[pos[0], pos[1]] = -spins[pos[0], pos[1]]
    return spins

#Función que genera los puntos para ser usado en la animación de ising
@jit(nopython=True)
def ising (N, T, n_data, dt):
    step = dt * N**2

    spins = init_spins(N)
    all_spins = np.zeros((n_data, N, N))

    i =0
    while (i < n_data):
        j = 0
        mag = 0
        while (j < step):
            spins = moves_spin(spins, N, T)
            j+=1
        all_spins[i] = spins
        i+=1
    return all_spins

i = 0
while (i < len(T)):
    ising_data[i] = ising(N, T[i], n_data, dt_ising)
    i +=1

def update(j_frame, frames_data, im):
    # Actualiza el gráfico con la configuración del sistema
    im.set_data(frames_data[j_frame])
    return im,

fig=plt.figure(figsize=(7,7)) #Size of the plot
ax=fig.add_subplot(111)

FFwriter = FFMpegWriter(fps=60, extra_args=['-vcodec', 'libx264'])

i = 0
while (i < len(T)):
    im = ax.imshow(ising_data[i][0], cmap="binary", vmin=-1, vmax=+1)
    animation = FuncAnimation(
            fig, update,
            fargs=(ising_data[i], im), frames=len(ising_data[i]), blit=True, interval=10)
    animation.save("{}.mp4".format(f"obligatoria2/test_time/video/ffm_N_{N}T_{(i+1)*0.5}"), dpi=100, writer = FFwriter)
    i+=1

plt.close()

stop = timeit.default_timer()
timer = stop - start

f = open("obligatoria2/test_time/runtime.txt", "a")
f.write(f"{N}, {timer}\n")
f.close()
