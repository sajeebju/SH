//This is the extended version code of - 
// https://github.com/daniel-koehn/Theory-of-seismic-waves-II/blob/master/06_2D_SH_Love_wave_modelling/2_From_2D_acoustic_to_SH_FD_modelling_final.ipynb

#include<iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <string>
#include <utility> 

class VelocityUpdater {
public:
    VelocityUpdater(double dx, double dz, double dt, int nx, int nz, const std::vector<std::vector<double>>& rho)
        : dx(dx), dz(dz), dt(dt), nx(nx), nz(nz), rho(rho) {}

    void update_velocity(std::vector<std::vector<double>>& vy,
                          const std::vector<std::vector<double>>& syx,
                          const std::vector<std::vector<double>>& syz,
                          int order) {
        switch (order) {
            case 2:
                update_vel_2nd_order(vy, syx, syz);
                break;
            case 4:
                update_vel_4th_order(vy, syx, syz);
                break;
            case 6:
                update_vel_6th_order(vy, syx, syz);
                break;
            case 8:
                update_vel_8th_order(vy, syx, syz);
                break;
            default:
                throw std::invalid_argument("Unsupported finite-difference order. Choose 2, 4, 6, or 8.");
        }
    }

private:
    double dx, dz, dt;
    int nx, nz;
    std::vector<std::vector<double>> rho;
    void update_vel_2nd_order(std::vector<std::vector<double>>& vy,
                           const std::vector<std::vector<double>>& syx,
                           const std::vector<std::vector<double>>& syz) {
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < nz - 1; ++j) {
                double syx_x = (syx[i][j] - syx[i - 1][j]) / dx;
                double syz_z = (syz[i][j] - syz[i][j - 1]) / dz;
                vy[i][j] += (dt / rho[i][j]) * (syx_x + syz_z);
            }
        }
    }

    void update_vel_4th_order(std::vector<std::vector<double>>& vy,
                           const std::vector<std::vector<double>>& syx,
                           const std::vector<std::vector<double>>& syz) {
        for (int i = 2; i < nx - 2; ++i) {
            for (int j = 2; j < nz - 2; ++j) {
                double syx_x = (-syx[i - 2][j] + 8 * syx[i - 1][j] - 8 * syx[i + 1][j] + syx[i + 2][j]) / (12 * dx);
                double syz_z = (-syz[i][j - 2] + 8 * syz[i][j - 1] - 8 * syz[i][j + 1] + syz[i][j + 2]) / (12 * dz);
                vy[i][j] += (dt / rho[i][j]) * (syx_x + syz_z);
            }
        }
    }

    void update_vel_6th_order(std::vector<std::vector<double>>& vy,
                           const std::vector<std::vector<double>>& syx,
                           const std::vector<std::vector<double>>& syz) {
        for (int i = 3; i < nx - 3; ++i) {
            for (int j = 3; j < nz - 3; ++j) {
                double syx_x = (-syx[i - 3][j] + 9 * syx[i - 2][j] - 45 * syx[i - 1][j] +
                                 45 * syx[i + 1][j] - 9 * syx[i + 2][j] + syx[i + 3][j]) / (60 * dx);
                double syz_z = (-syz[i][j - 3] + 9 * syz[i][j - 2] - 45 * syz[i][j - 1] +
                                 45 * syz[i][j + 1] - 9 * syz[i][j + 2] + syz[i][j + 3]) / (60 * dz);
                vy[i][j] += (dt / rho[i][j]) * (syx_x + syz_z);
            }
        }
    }

    void update_vel_8th_order(std::vector<std::vector<double>>& vy,
                           const std::vector<std::vector<double>>& syx,
                           const std::vector<std::vector<double>>& syz) {
        for (int i = 4; i < nx - 4; ++i) {
            for (int j = 4; j < nz - 4; ++j) {
                double syx_x = (-syx[i - 4][j] + 8 * syx[i - 3][j] - 28 * syx[i - 2][j] +
                                 56 * syx[i - 1][j] - 56 * syx[i + 1][j] + 28 * syx[i + 2][j] -
                                 8 * syx[i + 3][j] + syx[i + 4][j]) / (280 * dx);
                double syz_z = (-syz[i][j - 4] + 8 * syz[i][j - 3] - 28 * syz[i][j - 2] +
                                 56 * syz[i][j - 1] - 56 * syz[i][j + 1] + 28 * syz[i][j + 2] -
                                 8 * syz[i][j + 3] + syz[i][j + 4]) / (280 * dz);
                vy[i][j] += (dt / rho[i][j]) * (syx_x + syz_z);
            }
        }
    }
};

class StressUpdater {
public:
 
    StressUpdater(double dx, double dz, double dt, int nx, int nz, 
                  const std::vector<std::vector<double>>& mux, 
                  const std::vector<std::vector<double>>& muz)
        : dx(dx), dz(dz), dt(dt), nx(nx), nz(nz), mux(mux), muz(muz) {}

    void update_stress(std::vector<std::vector<double>>& syx,
                       std::vector<std::vector<double>>& syz,
                       const std::vector<std::vector<double>>& vy,
                       int order) {
        switch (order) {
            case 2:
                update_stress_2nd_order(syx, syz, vy);
                break;
            case 4:
                update_stress_4th_order(syx, syz, vy);
                break;
            case 6:
                update_stress_6th_order(syx, syz, vy);
                break;
            case 8:
                update_stress_8th_order(syx, syz, vy);
                break;
            default:
                throw std::invalid_argument("Unsupported finite-difference order. Choose 2, 4, 6, or 8.");
        }
    }

private:
    double dx, dz, dt;
    int nx, nz;
    std::vector<std::vector<double>> mux, muz;

    void update_stress_2nd_order(std::vector<std::vector<double>>& syx,
                                  std::vector<std::vector<double>>& syz,
                                  const std::vector<std::vector<double>>& vy) {
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < nz - 1; ++j) {
                double vy_x = (vy[i + 1][j] - vy[i][j]) / dx;
                double vy_z = (vy[i][j + 1] - vy[i][j]) / dz;
                syx[i][j] += dt * mux[i][j] * vy_x;
                syz[i][j] += dt * muz[i][j] * vy_z;
            }
        }
    }

    void update_stress_4th_order(std::vector<std::vector<double>>& syx,
                                  std::vector<std::vector<double>>& syz,
                                  const std::vector<std::vector<double>>& vy) {
        for (int i = 2; i < nx - 2; ++i) {
            for (int j = 2; j < nz - 2; ++j) {
                double vy_x = (-vy[i - 2][j] + 8 * vy[i - 1][j] - 8 * vy[i + 1][j] + vy[i + 2][j]) / (12 * dx);
                double vy_z = (-vy[i][j - 2] + 8 * vy[i][j - 1] - 8 * vy[i][j + 1] + vy[i][j + 2]) / (12 * dz);
                syx[i][j] += dt * mux[i][j] * vy_x;
                syz[i][j] += dt * muz[i][j] * vy_z;
            }
        }
    }

    void update_stress_6th_order(std::vector<std::vector<double>>& syx,
                                  std::vector<std::vector<double>>& syz,
                                  const std::vector<std::vector<double>>& vy) {
        for (int i = 3; i < nx - 3; ++i) {
            for (int j = 3; j < nz - 3; ++j) {
                double vy_x = (-vy[i - 3][j] + 9 * vy[i - 2][j] - 45 * vy[i - 1][j] +
                                45 * vy[i + 1][j] - 9 * vy[i + 2][j] + vy[i + 3][j]) / (60 * dx);
                double vy_z = (-vy[i][j - 3] + 9 * vy[i][j - 2] - 45 * vy[i][j - 1] +
                                45 * vy[i][j + 1] - 9 * vy[i][j + 2] + vy[i][j + 3]) / (60 * dz);
                syx[i][j] += dt * mux[i][j] * vy_x;
                syz[i][j] += dt * muz[i][j] * vy_z;
            }
        }
    }

   
    void update_stress_8th_order(std::vector<std::vector<double>>& syx,
                                  std::vector<std::vector<double>>& syz,
                                  const std::vector<std::vector<double>>& vy) {
        for (int i = 4; i < nx - 4; ++i) {
            for (int j = 4; j < nz - 4; ++j) {
                double vy_x = (-vy[i - 4][j] + 8 * vy[i - 3][j] - 28 * vy[i - 2][j] + 56 * vy[i - 1][j] -
                                56 * vy[i + 1][j] + 28 * vy[i + 2][j] - 8 * vy[i + 3][j] + vy[i + 4][j]) / (280 * dx);
                double vy_z = (-vy[i][j - 4] + 8 * vy[i][j - 3] - 28 * vy[i][j - 2] + 56 * vy[i][j - 1] -
                                56 * vy[i][j + 1] + 28 * vy[i][j + 2] - 8 * vy[i][j + 3] + vy[i][j + 4]) / (280 * dz);
                syx[i][j] += dt * mux[i][j] * vy_x;
                syz[i][j] += dt * muz[i][j] * vy_z;
            }
        }
    }
};

class ShearAverage {
public:

    ShearAverage(int nx, int nz)
        : nx(nx), nz(nz) {}
    void compute_average(std::vector<std::vector<double>>& mux,
                         std::vector<std::vector<double>>& muz,
                         const std::vector<std::vector<double>>& mu,
                         int order) {
        switch (order) {
            case 2:
                shear_avg_2nd_order(mux, muz, mu);
                break;
            case 4:
                shear_avg_4th_order(mux, muz, mu);
                break;
            case 6:
                shear_avg_6th_order(mux, muz, mu);
                break;
            case 8:
                shear_avg_8th_order(mux, muz, mu);
                break;
            default:
                throw std::invalid_argument("Unsupported finite-difference order. Choose 2, 4, 6, or 8.");
        }
    }
    
private:
    int nx, nz;

    void shear_avg_2nd_order(std::vector<std::vector<double>>& mux,
                              std::vector<std::vector<double>>& muz,
                              const std::vector<std::vector<double>>& mu) {
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < nz - 1; ++j) {
                mux[i][j] = 2.0 / (1.0 / mu[i + 1][j] + 1.0 / mu[i][j]);
                muz[i][j] = 2.0 / (1.0 / mu[i][j + 1] + 1.0 / mu[i][j]);
            }
        }
    }

    void shear_avg_4th_order(std::vector<std::vector<double>>& mux,
                              std::vector<std::vector<double>>& muz,
                              const std::vector<std::vector<double>>& mu) {
        for (int i = 2; i < nx - 2; ++i) {
            for (int j = 2; j < nz - 2; ++j) {
                mux[i][j] = 2.0 / (1.0 / mu[i + 2][j] + 1.0 / mu[i + 1][j] + 1.0 / mu[i][j]);
                muz[i][j] = 2.0 / (1.0 / mu[i][j + 2] + 1.0 / mu[i][j + 1] + 1.0 / mu[i][j]);
            }
        }
    }

    void shear_avg_6th_order(std::vector<std::vector<double>>& mux,
                              std::vector<std::vector<double>>& muz,
                              const std::vector<std::vector<double>>& mu) {
        for (int i = 3; i < nx - 3; ++i) {
            for (int j = 3; j < nz - 3; ++j) {
                mux[i][j] = 2.0 / (1.0 / mu[i + 3][j] + 1.0 / mu[i + 2][j] +
                                   1.0 / mu[i + 1][j] + 1.0 / mu[i][j]);
                muz[i][j] = 2.0 / (1.0 / mu[i][j + 3] + 1.0 / mu[i][j + 2] +
                                   1.0 / mu[i][j + 1] + 1.0 / mu[i][j]);
            }
        }
    }

    void shear_avg_8th_order(std::vector<std::vector<double>>& mux,
                              std::vector<std::vector<double>>& muz,
                              const std::vector<std::vector<double>>& mu) {
        for (int i = 4; i < nx - 4; ++i) {
            for (int j = 4; j < nz - 4; ++j) {
                mux[i][j] = 2.0 / (1.0 / mu[i + 4][j] + 1.0 / mu[i + 3][j] +
                                   1.0 / mu[i + 2][j] + 1.0 / mu[i + 1][j] + 1.0 / mu[i][j]);
                muz[i][j] = 2.0 / (1.0 / mu[i][j + 4] + 1.0 / mu[i][j + 3] +
                                   1.0 / mu[i][j + 2] + 1.0 / mu[i][j + 1] + 1.0 / mu[i][j]);
            }
        }
    }
};

std::vector<std::vector<double>> absorb(int nx, int nz, int w, double a) {
    std::vector<double> coeff(w);
    for (int i = 0; i < w; ++i) {
        coeff[i] = std::exp(-(a * a * (w - i) * (w - i)));
    }
    std::vector<std::vector<double>> absorb_coeff(nx, std::vector<double>(nz, 1.0));
    for (int i = 0; i < w; ++i) {
        int ze = nz - i - 1;
        for (int j = 0; j < ze; ++j) {
            absorb_coeff[i][j] = coeff[i];
        }
    }
    for (int i = 0; i < w; ++i) {
        int ii = nx - i - 1;
        int ze = nz - i - 1;
        for (int j = 0; j < ze; ++j) {
            absorb_coeff[ii][j] = coeff[i];
        }
    }
    for (int j = 0; j < w; ++j) {
        int jj = nz - j - 1;
        int xb = j;
        int xe = nx - j;
        for (int i = xb; i < xe; ++i) {
            absorb_coeff[i][jj] = coeff[j];
        }
    }

    return absorb_coeff;
}

std::pair<std::vector<double>, std::vector<double>> SH_SEIS(int nx, int nz, int nt, double dx, double dz, double dt,
                                                           double f0, double t0, int isrc, int jsrc, int ir, int jr,
                                                           const std::vector<std::vector<double>>& vs,
                                                           const std::vector<std::vector<double>>& rho, int w, double a,
                                                           int accuracy) {
    if (accuracy != 2 && accuracy != 4 && accuracy != 6 && accuracy != 8) {
        throw std::invalid_argument("Unsupported accuracy level. Choose among 2, 4, 6, or 8.");
    }

    // Gaussian STF
    std::vector<double> src(nt);
    std::vector<double> time(nt);
    for (int it = 0; it < nt; ++it) {
        time[it] = it * dt;
        src[it] = -2.0 * (time[it] - t0) * (f0 * f0) * std::exp(-(f0 * f0) * (time[it] - t0) * (time[it] - t0));
    }
    std::vector<std::vector<double>> vy(nx, std::vector<double>(nz, 0.0));
    std::vector<std::vector<double>> syx(nx, std::vector<double>(nz, 0.0));
    std::vector<std::vector<double>> syz(nx, std::vector<double>(nz, 0.0));
    std::vector<std::vector<double>> mu(nx, std::vector<double>(nz, 0.0));
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < nz; ++j) {
            mu[i][j] = rho[i][j] * vs[i][j] * vs[i][j];
        }
    }
    std::vector<std::vector<double>> mux = mu;
    std::vector<std::vector<double>> muz = mu;
    ShearAverage update_shear_avg(nx, nz);
    if (accuracy == 2) {
        update_shear_avg.compute_average(mux, muz, mu, 2);
    } else if (accuracy == 4) {
        update_shear_avg.compute_average(mux, muz, mu, 4);
    } else if (accuracy == 6) {
        update_shear_avg.compute_average(mux, muz, mu, 6);
    } else if (accuracy == 8) {
        update_shear_avg.compute_average(mux, muz, mu, 8);
    }

    std::vector<double> seis(nt, 0.0);
    auto absorb_coeff = absorb(nx, nz, w, a);

    VelocityUpdater updater_vel(dx, dz, dt, nx, nz, rho);
    StressUpdater updater_stress(dx, dz, dt, nx, nz, mux, muz);

    for (int it = 0; it < nt; ++it) {
        if (accuracy == 2) {
              updater_vel.update_velocity(vy, syx, syz, 2);
            updater_stress.update_stress(syx, syz, vy, 2);
        } else if (accuracy == 4) {
              updater_vel.update_velocity(vy, syx, syz, 4);
            updater_stress.update_stress(syx, syz, vy, 4);
        } else if (accuracy == 6) {
              updater_vel.update_velocity(vy, syx, syz, 6);
            updater_stress.update_stress(syx, syz, vy, 6);
        } else if (accuracy == 8) {
             updater_vel.update_velocity(vy, syx, syx, 8);
            updater_stress.update_stress(syx, syz, vy, 8);
        }
      
        vy[isrc][jsrc] += dt * src[it] / (rho[isrc][jsrc] * dx * dz);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nz; ++j) {
                vy[i][j] *= absorb_coeff[i][j];
                syx[i][j] *= absorb_coeff[i][j];
                syz[i][j] *= absorb_coeff[i][j];
            }
        }

        seis[it] = vy[ir][jr];
    }
    return {time, seis};
}
void writeSeisToFile(const std::string& filename, const std::pair<std::vector<double>, std::vector<double>>& data) {
    const auto& time = data.first;
    const auto& seis = data.second;
    if (time.size() != seis.size()) {
        throw std::runtime_error("Time and seis vectors must have the same size.");
    }
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::ios_base::failure("Failed to open file: " + filename);
    }
    for (size_t i = 0; i < time.size(); ++i) {
        outFile << time[i] << " " << seis[i] << "\n";
    }
    outFile.close();
}

int main() {
 
    double dx = 1.0;  
    double dz = 1.0;  
    double dt = 0.001;

    double xmax = 500.0; 
    double zmax = 500.0; 
    double tmax = 0.502; 

    int nx = static_cast<int>(xmax / dx); 
    int nz = static_cast<int>(zmax / dz); 
    int nt = static_cast<int>(tmax / dt); 

    double vs0 = 580.0;  
    double rho0 = 1000.0; 

    std::vector<std::vector<double>> vs(nx, std::vector<double>(nz, vs0));   
    std::vector<std::vector<double>> rho(nx, std::vector<double>(nz, rho0)); 

    double xsrc = 250.0; 
    double zsrc = 250.0; 
    double xr = 330.0; 
    double zr = 330.0; 
    int isrc = static_cast<int>(xsrc / dx);
    int jsrc = static_cast<int>(zsrc / dz);
    int ir = static_cast<int>(xr / dx);
    int jr = static_cast<int>(zr / dz);
    double f0 = 40.0; 
    double t0 = 4.0/f0;  
    int w = 60; 
    double a = 0.0053;
    int accuracy = 6; 
    auto result = SH_SEIS(nx, nz, nt, dx, dz, dt, f0, t0, isrc, jsrc, ir, jr, vs, rho, w, a, accuracy);
    writeSeisToFile("seis.txt", result);

    return 0;
}
