# distutils: language = c++
# distutils: extra_compile_args = -std=c++11

from libcpp.vector cimport vector
from libcpp.utility cimport pair
import numpy as np
cimport numpy as np

cdef extern from "sh2dwave.h":
    pair[vector[vector[vector[double]]], vector[vector[double]]] wave_propagate(
        int nx, int nz, int nt, double dx, double dz, double dt,
        double f0, double t0, int xsrc, int zsrc,
        const vector[vector[double]]& vs,
        const vector[vector[double]]& rho,
        int accuracy)

def py_wave_propagate(int nx, int nz, int nt, double dx, double dz, double dt,
                      double f0, double t0, int xsrc, int zsrc,
                      np.ndarray[np.double_t, ndim=2] vs,
                      np.ndarray[np.double_t, ndim=2] rho, int accuracy):
    cdef vector[vector[double]] c_vs = <vector[vector[double]]>([list(row) for row in vs])
    cdef vector[vector[double]] c_rho = <vector[vector[double]]>([list(row) for row in rho])
    cdef pair[vector[vector[vector[double]]], vector[vector[double]]] result
    result = wave_propagate(nx, nz, nt, dx, dz, dt, f0, t0, xsrc, zsrc, c_vs, c_rho, accuracy)
    return (np.array(result.first, dtype=np.float64), np.array(result.second, dtype=np.float64))

