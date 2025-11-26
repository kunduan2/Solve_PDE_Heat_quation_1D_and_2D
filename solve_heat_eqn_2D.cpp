#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main() {
    // ------------------------------- 2D Heat Equation (Explicit FTCS) --------------------------

    // Physical & numerical parameters
    double alpha = 0.0010;      // thermal diffusivity
    double Lx = 1.0;            // domain length in x-direction
    int Nx = 100;               // number of grid points in x
    double dx = Lx / (Nx - 1); // spatial step in x

    double Ly = 1.0;            // domain length in y-direction
    int Ny = 100;               // number of grid points in y
    double dy = Ly / (Ny - 1); // spatial step in y

    double dt = 0.0001;        // time step
    double T  = 0.1;           // final simulation time
    int Nt = T / dt;            // number of time steps

    // Stability condition for explicit FTCS scheme
    double r = alpha * dt / (dx * dx);
    if (r > 0.5) {
        cout << "WARNING: Unstable scheme, r = " << r << endl;
    }

    // 2D arrays to store the current and updated temperature fields
    vector<vector<double>> u(Nx, vector<double>(Ny));
    vector<vector<double>> u_new(Nx, vector<double>(Ny));

    // Initial condition: u(x,y,0) = sin(3x) * sin(3y)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u[i][j] = sin(3.0 * (i * dx)) * sin(3.0 * (j * dy));
        }
    }

    // ----------------------- Time-stepping loop -----------------------
    for (int t = 0; t < Nt; t++) {
        // Update interior points using 5-point Laplacian
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                u_new[i][j] =
                    u[i][j] + r * (
                        u[i+1][j] +  // right neighbor
                        u[i-1][j] +  // left neighbor
                        u[i][j+1] +  // top neighbor
                        u[i][j-1] -  // bottom neighbor
                        4.0 * u[i][j] // central point
                    );
            }
        }

        // Apply Dirichlet boundary conditions (u = 0 on boundaries)
        for (int i = 0; i < Nx; i++) u_new[i][0] = u_new[i][Ny-1] = 0.0; 
        for (int j = 0; j < Ny; j++) u_new[0][j] = u_new[Nx-1][j] = 0.0;

        // Swap new field into current field
        u = u_new;
    }

    // ----------------------- Save final temperature field -----------------------
    ofstream fout("heat2d.dat");
    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = j * dy;
            fout << x << " " << y << " " << u[i][j] << "\n";
        }
        fout << "\n"; // blank line for gnuplot pm3d
    }
    fout.close();

    // Quick visualization using gnuplot (3D surface)
    system("gnuplot -e \"set pm3d; splot 'heat2d.dat' with pm3d; pause -1\"");

    return 0;
}
