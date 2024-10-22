#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <iomanip> // for setting precision

using namespace std;

// Function to print hub and authority values with the required format
void printValues(const vector<double>& authorities, const vector<double>& hubs, int iteration) {
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
    for (double& v : values) {
        v /= sum;
    }
}

void hitsAlgorithm(const vector<vector<int>>& outEdges, const vector<vector<int>>& inEdges, int iterations, double errorRate, vector<double>& hubs, vector<double>& authorities) {
    int N = outEdges.size();
    vector<double> new_hubs(N, 0.0), new_authorities(N, 0.0);

    for (int iter = 1; iter <= iterations || (iterations == 0 && computeError(hubs, new_hubs) > errorRate); ++iter) {
        // Update authority scores: authorities[i] depends on the hub scores of pages linking to it
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

        // Check for convergence
        // if (computeError(hubs, new_hubs) < errorRate && computeError(authorities, new_authorities) < errorRate) {
        //     break;
        // }

        printValues(new_authorities, new_hubs, iter);

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
    vector<vector<int>> outEdges(n);  // Outgoing edges
    vector<vector<int>> inEdges(n);   // Incoming edges

    for (int i = 0; i < m; ++i) {
        int u, v;
        infile >> u >> v;
        outEdges[u].push_back(v);  // Edge from u -> v
        inEdges[v].push_back(u);   // Edge to v from u (reverse for authority calculation)
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

    printValues(authorities, hubs, 0);

    hitsAlgorithm(outEdges, inEdges, iterations, errorRate, hubs, authorities);

    return 0;
}
