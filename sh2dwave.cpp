#include "sh2dwave.h"
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iostream>


    VelocityUpdater::VelocityUpdater(double dx, double dz, double dt, int nx, int nz, const std::vector<std::vector<double>>& rho)
        : dx(dx), dz(dz), dt(dt), nx(nx), nz(nz), rho(rho) {}

    void  VelocityUpdater::update_velocity(std::vector<std::vector<double>>& vy,
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


    void  VelocityUpdater::update_vel_2nd_order(std::vector<std::vector<double>>& vy,
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

    void  VelocityUpdater::update_vel_4th_order(std::vector<std::vector<double>>& vy,
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

    void  VelocityUpdater::update_vel_6th_order(std::vector<std::vector<double>>& vy,
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

    void  VelocityUpdater::update_vel_8th_order(std::vector<std::vector<double>>& vy,
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

    StressUpdater::StressUpdater(double dx, double dz, double dt, int nx, int nz, 
                  const std::vector<std::vector<double>>& mux, 
                  const std::vector<std::vector<double>>& muz)
        : dx(dx), dz(dz), dt(dt), nx(nx), nz(nz), mux(mux), muz(muz) {}

    void  StressUpdater::update_stress(std::vector<std::vector<double>>& syx,
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

    void  StressUpdater::update_stress_2nd_order(std::vector<std::vector<double>>& syx,
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

    void  StressUpdater::update_stress_4th_order(std::vector<std::vector<double>>& syx,
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

    void  StressUpdater::update_stress_6th_order(std::vector<std::vector<double>>& syx,
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

   
    void  StressUpdater::update_stress_8th_order(std::vector<std::vector<double>>& syx,
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

    ShearAverage::ShearAverage(int nx, int nz)
        : nx(nx), nz(nz) {}
        
    void  ShearAverage::compute_average(std::vector<std::vector<double>>& mux,
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
    
    void  ShearAverage::shear_avg_2nd_order(std::vector<std::vector<double>>& mux,
                              std::vector<std::vector<double>>& muz,
                              const std::vector<std::vector<double>>& mu) {
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < nz - 1; ++j) {
                mux[i][j] = 2.0 / (1.0 / mu[i + 1][j] + 1.0 / mu[i][j]);
                muz[i][j] = 2.0 / (1.0 / mu[i][j + 1] + 1.0 / mu[i][j]);
            }
        }
    }

    void  ShearAverage::shear_avg_4th_order(std::vector<std::vector<double>>& mux,
                              std::vector<std::vector<double>>& muz,
                              const std::vector<std::vector<double>>& mu) {
        for (int i = 2; i < nx - 2; ++i) {
            for (int j = 2; j < nz - 2; ++j) {
                mux[i][j] = 2.0 / (1.0 / mu[i + 2][j] + 1.0 / mu[i + 1][j] + 1.0 / mu[i][j]);
                muz[i][j] = 2.0 / (1.0 / mu[i][j + 2] + 1.0 / mu[i][j + 1] + 1.0 / mu[i][j]);
            }
        }
    }

    void  ShearAverage::shear_avg_6th_order(std::vector<std::vector<double>>& mux,
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

    void  ShearAverage::shear_avg_8th_order(std::vector<std::vector<double>>& mux,
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




std::pair<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<double>>> wave_propagate(int nx, int nz, int nt, double dx, double dz, double dt,
                                                           double f0, double t0, int xsrc, int zsrc,
                                                           const std::vector<std::vector<double>>& vs,
                                                           const std::vector<std::vector<double>>& rho,
                                                           int accuracy) {


    if (accuracy != 2 && accuracy != 4 && accuracy != 6 && accuracy != 8) {
        throw std::invalid_argument("Unsupported accuracy level. Choose among 2, 4, 6, or 8.");
    }


        double isx = static_cast<int>(xsrc / dx);
        double isz = static_cast<int>(zsrc / dz);

    // Gaussian STF
    std::vector<double> ST(nt);
    for (int it = 0; it < nt; ++it) {
    double current_time = it * dt;
    ST[it] = -2.0 * (current_time - t0) * (f0 * f0) * std::exp(-(f0 * f0) * (current_time - t0) * (current_time - t0));
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

    std::vector<std::vector<std::vector<double>>> solution_space(nx, std::vector<std::vector<double>>(nz, std::vector<double>(nt, 0.0)));

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
             updater_vel.update_velocity(vy, syx, syz, 8);
            updater_stress.update_stress(syx, syz, vy, 8);
        }
    

   if (isx >= 0 && isx < nx && isz >= 0 && isz < nz) {
    vy[isx][isz] += dt * dt * ST[it] / rho[isx][isz];
}


        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nz; ++j) {
                solution_space[i][j][it] = vy[i][j];
            }
        }

    }

    return {solution_space, vy};
}











