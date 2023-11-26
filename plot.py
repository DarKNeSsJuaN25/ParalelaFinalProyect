
import matplotlib.pyplot as plt
import numpy as np

# Leer el camino desde el archivo
path = np.loadtxt("min_path.txt", dtype=int)

# Obtener las coordenadas de las ciudades desde el archivo
cities = np.loadtxt("output.txt", usecols=(1, 2))

# Graficar los puntos
plt.scatter(cities[:, 0], cities[:, 1], c='red', marker='o')

# Graficar el camino seguido por el viajero
for i in range(len(path) - 1):
    plt.plot([cities[path[i], 0], cities[path[i + 1], 0]],
             [cities[path[i], 1], cities[path[i + 1], 1]], c='blue')

# Conectar el último punto con el primero para cerrar el ciclo
plt.plot([cities[path[-1], 0], cities[path[0], 0]],
         [cities[path[-1], 1], cities[path[0], 1]], c='blue')

# Mostrar la gráfica
plt.savefig('grafica.png')
