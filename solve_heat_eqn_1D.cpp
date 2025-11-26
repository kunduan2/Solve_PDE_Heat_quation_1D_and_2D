#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main() {
    // ------------------------------- 1D Heat Equation (Explicit FTCS) --------------------------

    // Physical & numerical parameters
    double alpha = 0.0010;       // thermal diffusivity
    double L = 1.0;            // domain length
    int Nx = 1000;             // number of spatial grid points
    double dx = L / (Nx - 1);

    double dt = 0.0001;        // time step
    double T  = 0.1;           // final time
    int Nt = T / dt;

    // Stability condition
    double r = alpha * dt / (dx * dx);
    if (r > 0.5) {
        cout << "WARNING: Unstable scheme, r = " << r << endl;
    }

    vector<double> u(Nx), u_new(Nx);

    // Initial condition
    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        u[i] = sin(3.0 * x);
    }

    // Time stepping
    for (int t = 0; t < Nt; t++) {
        for (int i = 1; i < Nx - 1; i++) {
            u_new[i] = u[i] + r * (u[i+1] - 2*u[i] + u[i-1]);
        }

        // Dirichlet BCs
        u_new[0] = 0.0;
        u_new[Nx - 1] = 0.0;

        u = u_new;
    }

    // Save output
    ofstream fout("heat.dat");
    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        fout << x << " " << u[i] << "\n";
    }
    fout.close();

    // Plot using gnuplot
    system("gnuplot -e \"plot 'heat.dat' with lines lw 2; pause -1\"");

    return 0;
}
