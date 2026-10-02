// Verification test: pure Langmuir adsorption-desorption.
//
// With sticking coefficient 1 and constant desorption rate r_D, surface coverage obeys
//     d(theta)/dt = r_A (1 - theta) - r_D theta,
// whose solution from an empty surface is
//     theta(t) = r_A / (r_A + r_D) * (1 - exp(-(r_A + r_D) t)).
// On the lattice, r_A = (collisions per step) / (number of sites * dt).
//
// The test runs the stochastic solver, compares it with the analytical solution,
// writes both curves to langmuir_verification.txt, and fails if the error is too large.

#include <cmath>
#include <fstream>
#include <iostream>
#include "../src/solver.h"

int main() {
    solver_inputs inputs;
    inputs.langmuir_test = true;
    inputs.surface_size_X = 100;
    inputs.surface_size_Y = 100;
    inputs.time_step_size = 1e-3;   // s
    inputs.time_step_no = 2500;     // 2.5 s
    inputs.particle_flux = 16;      // gives r_A = 16 / (1e4 * 1e-3) = 1.6 1/s
    inputs.langmuir_k_des = 2.0;    // r_D = 2 1/s

    const double nSites = static_cast<double>(inputs.surface_size_X) * inputs.surface_size_Y;
    const double rA = inputs.particle_flux / (nSites * inputs.time_step_size);
    const double rD = inputs.langmuir_k_des;

    Solver solver(inputs);
    solver.execute();

    std::ofstream out("langmuir_verification.txt");
    out << "# time_s  theta_simulation  theta_analytical\n";

    double maxErr = 0.0;
    double sumSqErr = 0.0;
    const size_t n = solver.surface_cov.size();
    for (size_t i = 0; i < n; i++) {
        const double t = i * inputs.time_step_size;
        const double exact = rA / (rA + rD) * (1.0 - std::exp(-(rA + rD) * t));
        const double err = std::fabs(solver.surface_cov[i] - exact);
        maxErr = std::max(maxErr, err);
        sumSqErr += err * err;
        out << t << " " << solver.surface_cov[i] << " " << exact << "\n";
    }
    const double rmsErr = std::sqrt(sumSqErr / n);
    const double thetaSS = rA / (rA + rD);

    std::cout << "r_A = " << rA << " 1/s, r_D = " << rD << " 1/s\n"
              << "steady-state coverage: analytical " << thetaSS
              << ", simulation " << solver.surface_cov.back() << "\n"
              << "RMS error " << rmsErr << ", max error " << maxErr << "\n";

    // Statistical noise for 1e4 sites is about sqrt(theta(1-theta)/N) ~ 0.005
    const bool pass = rmsErr < 0.01 && maxErr < 0.03;
    std::cout << (pass ? "PASSED" : "FAILED") << std::endl;
    return pass ? 0 : 1;
}
