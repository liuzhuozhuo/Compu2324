# ================================================================================
# ANIMACION SISTEMA SOLAR 
#
# Genera una animación a partir de un fichero de datos con las posiciones
# de los planetas en diferentes instantes de tiempo.
# 
# El fichero debe estructurarse de la siguiente forma:
# 
#   x1_1, y1_1
#   x2_1, y2_1
#   x3_1, y3_1
#   (...)
#   xN_1, yN_1
#   
#   x1_2, y1_2
#   x2_2, y2_2
#   x3_2, y3_2
#   (...)
#   xN_2, yN_2
#
#   x1_3, y1_3
#   x2_3, y2_3
#   x3_3, y3_3
#   (...)
#   xN_3, yN_3
#   
#   (...)
#
# donde xi_j es la componente x del planeta i-ésimo en el instante de
# tiempo j-ésimo, e yi_j lo mismo en la componente y. El programa asume que
# el nº de planetas es siempre el mismo.
# ¡OJO! Los datos están separados por comas.
# 
# Si solo se especifica un instante de tiempo, se genera una imagen en pdf
# en lugar de una animación
#
# Se puede configurar la animación cambiando el valor de las variables
# de la sección "Parámetros"
#
# ================================================================================

# Importa los módulos necesarios
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Parámetros
# ========================================
file_in = f"obligatoria1/data/energy.txt" # Nombre del fichero de datos
file_out = f"obligatoria1/video/energy" # Nombre del fichero de salida (sin extensión)

planet_name = ['Sun', 'Mercury', 'Venus', 'Earth', 'Marth', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']


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

# El número de planetas es el número de líneas en cada bloque
# Lo calculamos del primer bloque
nplanets = len(frames_data[0])
frames_data = np.transpose(frames_data)


# Creación de la animación/gráfico
# ========================================
# Crea los objetos figure y axis

for i in range(nplanets):
    fig = plt.figure(figsize=(10, 5))    
    ax = fig.add_subplot(111) 
    plt.plot(frames_data[0][i], frames_data[1][i], label = "Cinética")
    plt.plot(frames_data[0][i], frames_data[2][i], label = "Potencial")
    plt.plot(frames_data[0][i], frames_data[3][i], label = "Cinética + Potencial 1")
    plt.plot(frames_data[0][i], frames_data[1][i] + frames_data[2][i], label = "Cinética + Potencial 2")
    ax.set_title(f"Energía para {planet_name[i]} (tiempo total: 2 años)")
    ax.grid()
    plt.legend()
    fig.savefig("{}.png".format(f'obligatoria1/plots/energy_{planet_name[i]}'))

data = pd.read_csv('obligatoria1/data/energy_total.txt', delimiter=",", header=0, names=["t", "K", "V", "TV"])

fig=plt.figure(figsize=(10,6)) #Size of the plot
ax=fig.add_subplot(111)
plt.plot(data.t, data.K, 'r-', label='Kinetic',linewidth=2)
plt.plot(data.t, data.V, 'b-', label='Potential',linewidth=2)
plt.plot(data.t, data.TV, 'g-', label='Kinetic + Potential',linewidth=2)
plt.xlabel('t ()',fontsize=12)
plt.ylabel('E ()',fontsize=12)
plt.legend(loc='best',fontsize=12)
plt.grid()
fig.savefig("{}.png".format('obligatoria1/plots/energy_total'))