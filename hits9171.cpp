/* First Name: Suraj Kumar
Last Name: Ojha
LAST FOUR NJIT ID: 9171 */

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <iomanip>

using namespace std;

// Function to print hub and authority values
void printValues(const vector<double>& authorities, const vector<double>& hubs, int iteration, bool largeGraph) {
    if (!largeGraph) {
        if (iteration == 0) {
            cout << "Base : " << iteration << " :";
        } else {
            cout << "Iter : " << iteration << " :";
        }
        for (size_t i = 0; i < authorities.size(); ++i) {
            cout << fixed << setprecision(7);
            cout << " A/H[ " << i << "]=" << authorities[i] << "/" << hubs[i];
        }
        cout << endl;
    } else if (largeGraph && iteration > 0) {
        cout << "Iter   : " << iteration << endl;
        for (size_t i = 0; i < authorities.size(); ++i) {
            cout << fixed << setprecision(7);
            cout << "A/H[ " << i << "]=" << authorities[i] << "/" << hubs[i] << endl;
        }
    }
}

double computeError(const vector<double>& old_vals, const vector<double>& new_vals) {
    double error = 0.0;
    for (size_t i = 0; i < old_vals.size(); ++i) {
        error = max(error, fabs(old_vals[i] - new_vals[i]));
    }
    return error;
}

void normalize(vector<double>& values) {
    double sum = 0.0;
    for (double v : values) {
        sum += v * v;
    }
    sum = sqrt(sum);

    // Avoid division by zero
    if (sum > 0) {
        for (double& v : values) {
            v /= sum;
        }
    }
    // If sum is zero, all values remain zero, avoiding NaN
}

void hitsAlgorithm(const vector<vector<int>>& outEdges, const vector<vector<int>>& inEdges, int iterations, double errorRate, vector<double>& hubs, vector<double>& authorities, bool largeGraph) {
    int N = outEdges.size();
    vector<double> new_hubs(N, 0.0), new_authorities(N, 0.0);

    if (!largeGraph) {
        printValues(authorities, hubs, 0, largeGraph);
    }

    int iter = 1;
    while (true) {
        for (int i = 0; i < N; ++i) {
            new_authorities[i] = 0.0;
            for (int neighbor : inEdges[i]) {
                new_authorities[i] += hubs[neighbor];
            }
        }

        // Update hub scores: hubs[i] depends on the authority scores of pages it links to
        for (int i = 0; i < N; ++i) {
            new_hubs[i] = 0.0;
            for (int neighbor : outEdges[i]) {
                new_hubs[i] += authorities[neighbor];
            }
        }

        // Normalize the hub and authority values
        normalize(new_hubs);
        normalize(new_authorities);

        if (!largeGraph) {
            printValues(new_authorities, new_hubs, iter, largeGraph);
        }

        // Check for errorRate
        if (iterations <= 0) {
            if (computeError(hubs, new_hubs) < errorRate && computeError(authorities, new_authorities) < errorRate) {
                break;
            }
        } else if (iter >= iterations) {
            break;
        }

        // Update the values for the next iteration
        hubs = new_hubs;
        authorities = new_authorities;
        ++iter;
    }

    // Print final values for large graphs
    if (largeGraph) {
        printValues(new_authorities, new_hubs, iter, largeGraph);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Please input the iteration/errorRate, initial values and graph file name\n";
        return 1;
    }

    int iterations = stoi(argv[1]);
    double errorRate = (iterations < 0) ? pow(10, iterations) : numeric_limits<double>::max();
    int initialization = stoi(argv[2]);
    string filename = argv[3];

    // read number of vertices and edges from file
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening file: " << filename << "\n";
        return 1;
    }

    int n, m;
    infile >> n >> m;

    // Determine if the graph is large
    bool largeGraph = (n > 10);

    // Adjust parameters for large graphs
    if (largeGraph) {
        errorRate = pow(10, -5);  // Fixed error rate
        initialization = -1;     // Initial values 1/n
        iterations = -1;         // Ignore user-provided iterations
    }

    // Create two separate vectors for incoming and outgoing edges
    vector<vector<int>> outEdges(n);
    vector<vector<int>> inEdges(n);

    for (int i = 0; i < m; ++i) {
        int u, v;
        infile >> u >> v;
        outEdges[u].push_back(v);  // For outgoing
        inEdges[v].push_back(u);   // For incoming
    }

    // Initialize hub and authority values
    vector<double> hubs(n, 0.0), authorities(n, 0.0);
    if (initialization == 1) {
        fill(hubs.begin(), hubs.end(), 1.0);
        fill(authorities.begin(), authorities.end(), 1.0);
    } else if (initialization == -1) {
        double init_val = 1.0 / n;
        fill(hubs.begin(), hubs.end(), init_val);
        fill(authorities.begin(), authorities.end(), init_val);
    } else if (initialization == -2) {
        double init_val = 1.0 / sqrt(n);
        fill(hubs.begin(), hubs.end(), init_val);
        fill(authorities.begin(), authorities.end(), init_val);
    }

    hitsAlgorithm(outEdges, inEdges, iterations, errorRate, hubs, authorities, largeGraph);

    // Close the file
    infile.close();

    return 0;
}
