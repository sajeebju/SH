#include "sh2dwave.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>

void writeSeisToFile(const std::string& filename, const std::pair<std::vector<double>, std::vector<double>>& data) {
    const auto& time = data.first;
    const auto& seis = data.second;
   
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::ios_base::failure("Failed to open file: " + filename);
    }
    for (size_t i = 0; i < time.size(); ++i) {
        outFile << time[i] << " " << seis[i] << "\n";
    }
    outFile.close();
    std::cout << "Data successfully written to " << filename << std::endl; // Debug statement
}

int main() {

    double dx = 1.0, dz = 1.0, dt = 0.001;
    double xmax = 500.0, zmax = 500.0, tmax = 0.502;
    int nx = static_cast<int>(xmax / dx);
    int nz = static_cast<int>(zmax / dz);
    int nt = static_cast<int>(tmax / dt);
    double vs0 = 580.0, rho0 = 1000.0;

    std::vector<std::vector<double>> vs(nx, std::vector<double>(nz, vs0));
    std::vector<std::vector<double>> rho(nx, std::vector<double>(nz, rho0));

    double xsrc = 250.0, zsrc = 250.0, xr = 330.0, zr = 330.0;
    int isrc = static_cast<int>(xsrc / dx);
    int jsrc = static_cast<int>(zsrc / dz);
    int ir = static_cast<int>(xr / dx);
    int jr = static_cast<int>(zr / dz);
    double f0 = 40.0, t0 = 4.0 / f0;
    int w = 60;
    double a = 0.0053;
    int accuracy = 6;

    auto result = SH_SEIS(nx, nz, nt, dx, dz, dt, f0, t0, isrc, jsrc, ir, jr, vs, rho, w, a, accuracy);
    writeSeisToFile("seismogram.txt", result);
    std::cout << "Seismogram calculated successfully." << std::endl; 

    return 0;
}

