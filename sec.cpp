#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <chrono>
#include <queue>
#include <omp.h>
#include <iomanip>
#include <malloc.h>

long int getMemoryUsage() {
    struct mallinfo2 info = mallinfo2();
    return info.uordblks / (1024 * 1024); // Total size of allocated blocks
}

struct City {
    int id;
    double x, y;
};

struct Node {
    int vertex;
    std::vector<std::pair<int, int>> path;
    std::vector<std::vector<int>> matrix_reduced;
    double cost;
    int level;
};

const int INF = std::numeric_limits<int>::max();

std::vector<std::vector<int>> createGraph(const std::vector<City>& cities) {
    int n = cities.size();
    std::vector<std::vector<int>> graph(n, std::vector<int>(n, INF));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double distance = std::sqrt(std::pow(cities[i].x - cities[j].x, 2) + std::pow(cities[i].y - cities[j].y, 2));
                graph[i][j] = static_cast<int>(distance);
            }
        }
    }

    return graph;
}

class CompareNodes {
public:
    bool operator()(const Node* lhs, const Node* rhs) const {
        return lhs->cost > rhs->cost;
    }
};

Node* newNode(const std::vector<std::vector<int>>& matrix_parent, const std::vector<std::pair<int, int>>& path, int level, int i, int j, int N) {
    auto node = new Node;
    node->path = path;
    if (level != 0)
        node->path.push_back(std::make_pair(i, j));
    node->matrix_reduced = matrix_parent;
    for (int k = 0; level != 0 && k < N; k++) {
        node->matrix_reduced[i][k] = INF;
        node->matrix_reduced[k][j] = INF;
    }

    node->matrix_reduced[j][0] = INF;
    node->level = level;
    node->vertex = j;
    return node;
}

void reduceRow(std::vector<std::vector<int>>& matrix_reduced, std::vector<int>& row, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
            if (matrix_reduced[i][j] < row[i])
                row[i] = matrix_reduced[i][j];

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
            if (matrix_reduced[i][j] != INF && row[i] != INF)
                matrix_reduced[i][j] -= row[i];
}

void reduceColumn(std::vector<std::vector<int>>& matrix_reduced, std::vector<int>& col, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
            if (matrix_reduced[i][j] < col[j])
                col[j] = matrix_reduced[i][j];

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
            if (matrix_reduced[i][j] != INF && col[j] != INF)
                matrix_reduced[i][j] -= col[j];
}

int costCalculation(std::vector<std::vector<int>>& matrix_reduced, int N) {
    int cost = 0;
    std::vector<int> row(N, INF);
    reduceRow(matrix_reduced, row, N);
    std::vector<int> col(N, INF);
    reduceColumn(matrix_reduced, col, N);

    #pragma omp parallel for reduction(+:cost)
    for (int i = 0; i < N; i++)
        cost += (row[i] != INF) ? row[i] : 0, cost += (col[i] != INF) ? col[i] : 0;

    return cost;
}

void printPath(const std::vector<std::pair<int, int>>& list) {
    for (int i = 0; i < list.size(); i++)
        std::cout << list[i].first + 1 << " -> " << list[i].second + 1 << std::endl;
}

Node* TSPBranchAndBound(const std::vector<std::vector<int>>& adjacencyMatrix) {
    std::priority_queue<Node*, std::vector<Node*>, CompareNodes> pq;
    std::vector<std::pair<int, int>> v;
    int N = adjacencyMatrix[0].size();

    auto root = newNode(adjacencyMatrix, v, 0, -1, 0, N);

    root->cost = costCalculation(root->matrix_reduced, N);

    pq.push(root);

    while (!pq.empty()) {
        auto min = pq.top();
        pq.pop();
        int i = min->vertex;
        if (min->level == N - 1) {
            min->path.push_back(std::make_pair(i, 0));
            return min;
        }

        #pragma omp parallel for shared(pq) schedule(dynamic)
        for (int j = 0; j < N; j++) {
            if (min->matrix_reduced[i][j] != INF) {
                Node* child = newNode(min->matrix_reduced, min->path, min->level + 1, i, j, N);
                child->cost = min->cost + min->matrix_reduced[i][j] + costCalculation(child->matrix_reduced, N);
                #pragma omp critical
                pq.push(child);
            }
        }

        delete min;
    }
    return nullptr;
}

void writeMinPathToFile(const std::vector<std::pair<int, int>>& path) {
    std::ofstream file("min_path.txt");
    for (const auto& pair : path) {
        file << pair.first << " " << pair.second << std::endl;
    }
    file.close();
}

int main() {
    omp_set_num_threads(4);
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

    std::vector<std::vector<int>> adjacencyMatrix = createGraph(cities);
    auto start = std::chrono::high_resolution_clock::now();
    auto ans = TSPBranchAndBound(adjacencyMatrix);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    std::cout << "Elapsed Time: " << std::fixed << std::setprecision(5) << elapsed.count() << " seconds" << std::endl;
    std::cout << "Memory Usage: " << getMemoryUsage() << " MB" << std::endl;

    // Camino minimo
    if (ans) {
        writeMinPathToFile(ans->path);
        std::cout << "Optimal Path:\n";
        printPath(ans->path);
        std::cout << "Cost: " << ans->cost << std::endl;
    } else {
        std::cout << "No solution found." << std::endl;
    }

    std::ofstream pythonFile("adjacency_matrix.txt");
    for (const auto& row : adjacencyMatrix) {
        for (int distance : row) {
            pythonFile << distance << " ";
        }
        pythonFile << std::endl;
    }
    pythonFile.close();

    // Crear un script de Python para graficar
    std::ofstream pythonScript("plot.py");
    pythonScript << R"(
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
)";
    pythonScript.close();

    // Ejecutar el script de Python
    system("python plot.py");

    return 0;
}
