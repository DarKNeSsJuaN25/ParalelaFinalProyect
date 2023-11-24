#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <chrono>
#include <omp.h>

struct City {
    int id;
    double x, y;
};

struct Node {
    int cityId;
    std::vector<int> path;
    double cost;
};

double euclideanDistance(const City& city1, const City& city2) {
    double dx = city1.x - city2.x;
    double dy = city1.y - city2.y;
    return std::sqrt(dx * dx + dy * dy);
}

double calculateCost(const std::vector<int>& path, const std::vector<City>& cities) {
    double cost = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        cost += euclideanDistance(cities[path[i]], cities[path[i + 1]]);
    }
    // Regresar al inicio
    cost += euclideanDistance(cities[path.back()], cities[path.front()]);
    return cost;
}

bool compareNodes(const Node& node1, const Node& node2) {
    return node1.cost > node2.cost;  // Ordenar en orden descendente
}

void branchAndBound(const std::vector<City>& cities, bool& finished, omp_lock_t& lock) {
    int n = cities.size();

    // Inicializar la raíz
    Node root;
    root.cityId = 0;  // Empezar desde la primera ciudad
    root.path.push_back(0);
    root.cost = std::numeric_limits<double>::infinity();  // Usar infinito como cota superior inicial

    // Crear la cola de prioridad para nodos
    std::vector<Node> priorityQueue;
    priorityQueue.push_back(root);

    // Variables para la mejor solución encontrada por cualquier hilo
    Node bestSolution;
    bestSolution.cost = std::numeric_limits<double>::infinity();

    // Variable para indicar que no hay más nodos que procesar
    finished = false;

    #pragma omp parallel
    {
        std::vector<Node> privateQueue;  // Cola privada para cada hilo

        while (true) {
            Node currentNode;

            #pragma omp critical
            {
                if (!priorityQueue.empty()) {
                    currentNode = priorityQueue.back();
                    priorityQueue.pop_back();
                } else {
                    if (!privateQueue.empty()) {
                        currentNode = privateQueue.back();
                        privateQueue.pop_back();
                    } else {
                        #pragma omp flush(finished)
                        if (finished) {
                            // No hay más nodos que procesar
                            break;
                        }
                    }
                }
            }

            // Resto del código para generar hijos y actualizar la mejor solución
            if (currentNode.path.size() == n) {
                double currentCost = calculateCost(currentNode.path, cities);

                // Actualizar la mejor solución del hilo si es necesario
                #pragma omp critical
                {
                    if (currentCost < bestSolution.cost) {
                        bestSolution = currentNode;
                    }
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    if (std::find(currentNode.path.begin(), currentNode.path.end(), i) == currentNode.path.end()) {
                        Node child;
                        child.cityId = i;
                        child.path = currentNode.path;
                        child.path.push_back(i);
                        child.cost = calculateCost(child.path, cities);

                        // Agregar a la cola privada del hilo
                        privateQueue.push_back(child);
                    }
                }
            }
        } // Fin del bucle en paralelo
    }

    // Mostrar la mejor solución
    std::cout << "Mejor recorrido: ";
    for (int cityId : bestSolution.path) {
        std::cout << cityId + 1 << " ";  // Sumar 1 para convertir de 0-indexed a 1-indexed
    }
    std::cout << bestSolution.path.front() + 1 <<  std::endl;
    std::cout << "Costo mínimo: " << calculateCost(bestSolution.path, cities) << std::endl;
}

int main() {
    std::ifstream file("output.txt");
    omp_set_num_threads(2);
    if (!file) {
        std::cerr << "Error al abrir el archivo." << std::endl;
        return 1;
    }

    std::vector<City> cities;
    int id, x, y;

    while (file >> id >> x >> y) {
        City city;
        city.id = id;
        city.x = x;
        city.y = y;
        cities.push_back(city);
    }

    file.close();

    // Inicializar el lock
    omp_lock_t lock;
    omp_init_lock(&lock);

    // Variable para indicar que no hay más nodos que procesar
    bool finished;

    auto start_time = std::chrono::high_resolution_clock::now();
    branchAndBound(cities, finished, lock);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    std::cout << "Tiempo de ejecución: " << duration.count() << " segundos" << std::endl;

    // Liberar el lock
    omp_destroy_lock(&lock);

    return 0;
}
