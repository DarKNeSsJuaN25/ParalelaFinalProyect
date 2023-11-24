#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <chrono>
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

void branchAndBound(const std::vector<City>& cities) {
    int n = cities.size();
    
    // Inicializar la raíz
    Node root;
    root.cityId = 0;  // Empezar desde la primera ciudad
    root.path.push_back(0);
    root.cost = std::numeric_limits<double>::infinity();  // Usar infinito como cota superior inicial

    // Crear la cola de prioridad para nodos
    std::vector<Node> priorityQueue;
    priorityQueue.push_back(root);

    while (!priorityQueue.empty()) {
        // Obtener el nodo con el menor costo
        std::pop_heap(priorityQueue.begin(), priorityQueue.end(), compareNodes);
        Node currentNode = priorityQueue.back();
        priorityQueue.pop_back();

        // Si el nodo es una hoja, actualizar la mejor solución si es necesario
        if (currentNode.path.size() == n) {
            double currentCost = calculateCost(currentNode.path, cities);
            if (currentCost < root.cost) {
                root = currentNode;
            }
        }

        // Generar hijos y agregarlos a la cola de prioridad
        for (int i = 0; i < n; ++i) {
            if (std::find(currentNode.path.begin(), currentNode.path.end(), i) == currentNode.path.end()) {
                Node child;
                child.cityId = i;
                child.path = currentNode.path;
                child.path.push_back(i);
                child.cost = calculateCost(child.path, cities);

                // Agregar a la cola de prioridad
                priorityQueue.push_back(child);
                std::push_heap(priorityQueue.begin(), priorityQueue.end(), compareNodes);
            }
        }
    }

    // Mostrar la mejor solución
    std::cout << "Mejor recorrido: ";
    for (int cityId : root.path) {
        std::cout << cityId + 1 << " ";  // Sumar 1 para convertir de 0-indexed a 1-indexed
    }
    std::cout << root.path.front() + 1 <<  std::endl;
    std::cout << "Costo mínimo: " << calculateCost(root.path, cities) << std::endl;
}

int main() {
    std::ifstream file("output.txt");

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
    // for(auto&i : cities){
    //     std::cout << i.id << " " << i.x << " "<< i.y << std::endl;
    // }
    auto start_time = std::chrono::high_resolution_clock::now();
    branchAndBound(cities);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    std::cout << "Tiempo de ejecución: " << duration.count() << " segundos" << std::endl;

    return 0;
}

