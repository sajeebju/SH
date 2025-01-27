// simulation.h
#ifndef SH2DWAVE_H
#define SH2DWAVE_H

#include <vector>
#include <string>
#include <utility>
#include<iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

class VelocityUpdater {
public:
    VelocityUpdater(double dx, double dz, double dt, int nx, int nz, const std::vector<std::vector<double>>& rho);
    void update_velocity(std::vector<std::vector<double>>& vy, const std::vector<std::vector<double>>& syx, const std::vector<std::vector<double>>& syz, int order);

private:
    double dx, dz, dt;
    int nx, nz;
    std::vector<std::vector<double>> rho;
    void update_vel_2nd_order(std::vector<std::vector<double>>& vy, const std::vector<std::vector<double>>& syx, const std::vector<std::vector<double>>& syz);
    void update_vel_4th_order(std::vector<std::vector<double>>& vy, const std::vector<std::vector<double>>& syx, const std::vector<std::vector<double>>& syz);
    void update_vel_6th_order(std::vector<std::vector<double>>& vy, const std::vector<std::vector<double>>& syx, const std::vector<std::vector<double>>& syz);
    void update_vel_8th_order(std::vector<std::vector<double>>& vy, const std::vector<std::vector<double>>& syx, const std::vector<std::vector<double>>& syz);
};

class StressUpdater {
public:
    StressUpdater(double dx, double dz, double dt, int nx, int nz, const std::vector<std::vector<double>>& mux, const std::vector<std::vector<double>>& muz);
    void update_stress(std::vector<std::vector<double>>& syx, std::vector<std::vector<double>>& syz, const std::vector<std::vector<double>>& vy, int order);

private:
    double dx, dz, dt;
    int nx, nz;
    std::vector<std::vector<double>> mux, muz;
    void update_stress_2nd_order(std::vector<std::vector<double>>& syx, std::vector<std::vector<double>>& syz, const std::vector<std::vector<double>>& vy);
    void update_stress_4th_order(std::vector<std::vector<double>>& syx, std::vector<std::vector<double>>& syz, const std::vector<std::vector<double>>& vy);
    void update_stress_6th_order(std::vector<std::vector<double>>& syx, std::vector<std::vector<double>>& syz, const std::vector<std::vector<double>>& vy);
    void update_stress_8th_order(std::vector<std::vector<double>>& syx, std::vector<std::vector<double>>& syz, const std::vector<std::vector<double>>& vy);
};

class ShearAverage {
public:
    ShearAverage(int nx, int nz);
    void compute_average(std::vector<std::vector<double>>& mux,
                         std::vector<std::vector<double>>& muz,
                         const std::vector<std::vector<double>>& mu,
                         int order);

private:
    int nx, nz;
    void shear_avg_2nd_order(std::vector<std::vector<double>>& mux,
                             std::vector<std::vector<double>>& muz,
                             const std::vector<std::vector<double>>& mu);
    void shear_avg_4th_order(std::vector<std::vector<double>>& mux,
                             std::vector<std::vector<double>>& muz,
                             const std::vector<std::vector<double>>& mu);
    void shear_avg_6th_order(std::vector<std::vector<double>>& mux,
                             std::vector<std::vector<double>>& muz,
                             const std::vector<std::vector<double>>& mu);
    void shear_avg_8th_order(std::vector<std::vector<double>>& mux,
                             std::vector<std::vector<double>>& muz,
                             const std::vector<std::vector<double>>& mu);
};

std::vector<std::vector<double>> absorb(int nx, int nz, int w, double a);

// Function to compute seismograms
std::pair<std::vector<double>, std::vector<double>> SH_SEIS(int nx, int nz, int nt, double dx, double dz, double dt,
                                                           double f0, double t0, int isrc, int jsrc, int ir, int jr,
                                                           const std::vector<std::vector<double>>& vs,
                                                           const std::vector<std::vector<double>>& rho, int w, double a,
                                                           int accuracy);
// Function to wave propagation
std::pair<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<double>>> wave_propagate(int nx, int nz, int nt, double dx, double dz, double dt,
                                                           double f0, double t0, int xsrc, int zsrc,
                                                           const std::vector<std::vector<double>>& vs,
                                                           const std::vector<std::vector<double>>& rho,
                                                           int accuracy);

#endif // SH2DWAVE_H
