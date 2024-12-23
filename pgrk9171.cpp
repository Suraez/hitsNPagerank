/* First Name: Suraj Kumar
Last Name: Ojha
LAST FOUR NJIT ID: 9171 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

void pageRankAlgorithm(int n, int m, vector<vector<int>>& adjList, int N, double computeError, int init) {
    vector<double> rank(n, 0.0);
    vector<double> newRank(n, 0.0);
    vector<int> outDegree(n, 0);
    double d = 0.85;  // Damping factor
    double initVal;

    bool largeGraph = (n > 10);

    // Initialize PageRank values
    switch (init) {
        case 0:
            initVal = 0.0;
            break;
        case 1:
            initVal = 1.0;
            break;
        case -1:
            initVal = 1.0 / n;
            break;
        case -2:
            initVal = 1.0 / sqrt(n);
            break;
        default:
            initVal = 0.0;
    }

    // Set initial PageRank values and compute out-degree for each node
    for (int i = 0; i < n; i++) {
        rank[i] = initVal;
        outDegree[i] = adjList[i].size();  // Out-degree
    }

    if (!largeGraph) {
        cout << "Base : 0 :";
        for (int i = 0; i < n; i++) {
            cout << " PR[ " << i << "]=" << fixed << setprecision(7) << rank[i];
        }
        cout << endl;
    }

    int iter = 0;
    bool stop = false;

    while (!stop) {
        iter++;

        for (int i = 0; i < n; i++) {
            newRank[i] = (1 - d) / n;
        }

        // Update PageRank values based on incoming links
        for (int i = 0; i < n; i++) {
            for (int j : adjList[i]) {
                if (outDegree[i] > 0) {
                    newRank[j] += d * rank[i] / outDegree[i];
                }
            }
        }

        // Check for convergence by calculating the error
        double maxError = 0.0;
        for (int i = 0; i < n; i++) {
            double error = fabs(newRank[i] - rank[i]);
            if (error > maxError) {
                maxError = error;
            }
        }

        if (!largeGraph) {
            cout << "Iter : " << iter << " :";
            for (int i = 0; i < n; i++) {
                cout << " PR[ " << i << "]=" << fixed << setprecision(7) << newRank[i];
            }
            cout << endl;
        }

        // Convergence conditions
        if (largeGraph) {
            if (maxError <= pow(10, -5)) {  // Fixed convergence threshold for large graphs
                stop = true;
            }
        } else {
            if (N > 0 && iter >= N) {
                stop = true;
            } else if (N == 0 && maxError <= pow(10, -5)) {
                stop = true;
            } else if (N < 0 && maxError <= pow(10, -abs(N))) {
                stop = true;
            }
        }

        rank = newRank;  // Update rank values for the next iteration
    }

    // Print final PageRank values if n > 10
    if (largeGraph) {
        cout << "Iter :" << iter << endl;
        for (int i = 0; i < n; i++) {
            cout << "P[" << i << "]=" << fixed << setprecision(7) << rank[i] << endl;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Please input the iteration/errorRate, initial values and graph file name\n";
        return 1;
    }

    int N = stoi(argv[1]);  // Number of iterations or error rate
    int init = stoi(argv[2]);  // Initialization values
    string filename = argv[3];  // filename

    ifstream file(filename);
    if (!file) {
        cout << "Error opening file: " << filename << endl;
        return 1;
    }

    int n, m;
    file >> n >> m;

    vector<vector<int>> adjList(n);  // Adjacency list

    int u, v;
    for (int i = 0; i < m; i++) {
        file >> u >> v;
        adjList[u].push_back(v);  // Directed edge from u to v
    }

    file.close();

    double computeError;
    if (n > 10) {
        computeError = pow(10, -5);  // Fixed error threshold for large graphs
        init = -1;  // Default initialization for large graphs
        N = -1;  // Ignore user-provided value of N
    } else {
        if (N == 0) {
            computeError = pow(10, -5);
        } else if (N < 0) {
            computeError = pow(10, -abs(N));
        } else {
            computeError = 0.0000001;
        }
    }

    pageRankAlgorithm(n, m, adjList, N, computeError, init);

    return 0;
}
