import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def read_snapshot(filepath):
    # Leer datos del archivo snapshot
    with open(filepath, 'r') as file:
        lines = file.readlines()
        # Extraer el tiempo del segundo dato de la primera línea que contiene metadatos
        time = float(lines[0].split()[1])
        # Saltar la primera línea que contiene metadatos
        lines = lines[1:]
        # Extraer las posiciones x e y y el tipo de cada partícula de cada línea
        data = np.array([[float(val) for val in line.split()[1:]] for line in lines])
    return time, data


def update(frame):
    ax.clear()
    ax.set_xlim(0, 1)  # Ajustar límites x según tus necesidades
    ax.set_ylim(0, 1)  # Ajustar límites y según tus necesidades
    ax.scatter(frame[1][:, 0], frame[1][:, 1], color='b')  # Dibujar las partículas en la instantánea actual
    # Seleccionar las partículas con itype = -1 y pintarlas de color gris
    gray_particles = frame[1][frame[1][:, 8] == -1]
    ax.scatter(gray_particles[:, 0], gray_particles[:, 1], color='gray')
    # Mostrar el tiempo en la gráfica
    ax.text(0.05, 0.95, f't = {frame[0]:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=12)


# Definir la carpeta donde se encuentran los archivos
folder = './'

# Crear una figura y un eje para la animación
fig, ax = plt.subplots()

# Lista para almacenar los datos de todas las instantáneas
data_frames = []

# Iterar sobre los nombres de archivo del tipo 'snapshot_0001', 'snapshot_0002', ..., 'snapshot_2000'
for i in range(1, 2001):
    # Formatear el nombre de archivo con el número de snapshot
    filename = f'snapshot_{i:04d}'
    # Crear la ruta completa al archivo
    filepath = os.path.join(folder, filename)
    # Verificar si el archivo existe
    if os.path.exists(filepath):
        # Leer datos del archivo snapshot y almacenarlos en la lista
        time, data = read_snapshot(filepath)
        data_frames.append((time, data))
    else:
        print(f'El archivo {filename} no existe en la carpeta especificada.')

# Crear la animación
anim = FuncAnimation(fig, update, frames=data_frames, interval=10)  # Intervalo de 50 milisegundos

# Guardar la animación como un archivo .mp4
anim.save('animacion02.mp4', writer='ffmpeg')

plt.show()
