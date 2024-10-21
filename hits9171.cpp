#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits>

using namespace std;

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
    for (double& v : values) {
        v /= sum;
    }
}

void hitsAlgorithm(const vector<vector<int>>& adjList, int iterations, double errorRate, vector<double>& hubs, vector<double>& authorities) {
    int N = adjList.size();
    vector<double> new_hubs(N, 0.0), new_authorities(N, 0.0);

    for (int iter = 0; iter < iterations || (iterations == 0 && computeError(hubs, new_hubs) > errorRate); ++iter) {
        // Update authority scores
        for (int i = 0; i < N; ++i) {
            new_authorities[i] = 0.0;
            for (int neighbor : adjList[i]) {
                new_authorities[i] += hubs[neighbor];
            }
        }

        // Update hub scores
        for (int i = 0; i < N; ++i) {
            new_hubs[i] = 0.0;
            for (int neighbor : adjList[i]) {
                new_hubs[i] += authorities[neighbor];
            }
        }

        // Normalize
        normalize(new_hubs);
        normalize(new_authorities);

        // Check for convergence
        if (computeError(hubs, new_hubs) < errorRate && computeError(authorities, new_authorities) < errorRate) {
            break;
        }

        hubs = new_hubs;
        authorities = new_authorities;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <iterations or errorrate> <initialization> <filename>\n";
        return 1;
    }

    int iterations = stoi(argv[1]);
    double errorRate = (iterations < 0) ? pow(10, iterations) : numeric_limits<double>::max();
    int initialization = stoi(argv[2]);
    string filename = argv[3];

    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening file: " << filename << "\n";
        return 1;
    }

    // Read the graph
    int n, m;
    infile >> n >> m;
    vector<vector<int>> adjList(n);

    for (int i = 0; i < m; ++i) {
        int u, v;
        infile >> u >> v;
        adjList[u].push_back(v);
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

    // Run the HITS algorithm
    hitsAlgorithm(adjList, iterations, errorRate, hubs, authorities);

    // Output the results
    cout << "Hub values:\n";
    for (double hub : hubs) {
        cout << hub << "\n";
    }
    cout << "Authority values:\n";
    for (double authority : authorities) {
        cout << authority << "\n";
    }

    return 0;
}
