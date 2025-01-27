#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from sh2dwave_wrapper import py_wave_propagate

xmax, zmax = 500.0, 500.0
nx, nz = 201, 201
Lx, Ly = xmax, zmax
dx, dz = Lx / (nx - 1), Ly / (nz - 1)
xsrc, zsrc = 250.0, 5.0 

rho = np.ones((nx, nz)) * 1000.0 
vs = np.ones((nx, nz)) * 580.0   

CN = 0.5  
tmax = 2.00  
nt = int(tmax / (CN * dx / np.max(vs)))
dt = CN * dx / np.max(vs)  
f0 = 40.0
t0 = 4.0 / 40.0  
accuracy = 2

solution_space, vy = py_wave_propagate(nx, nz, nt, dx, dz, dt, t0, f0, xsrc, zsrc, vs, rho, accuracy)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
k = 0
mskip = 5

vscale = (vy.T - np.mean(vy))
vscale = -vscale / np.max(vscale)

def SH_2D(i):
    global k
    wave_field = solution_space[:,:,k].T / np.max(np.abs(solution_space))
    data = wave_field + 0.0025 * vscale[:, ::-1]
    ax1.clear()
    plt.imshow(data, cmap='seismic', extent=[0, Lx, Ly, 0], vmin=-0.08, vmax=0.08)
    plt.xlim([0, Lx])
    plt.ylim([Ly, 0])
    plt.xlabel('Distance (m)', fontsize=12)
    plt.ylabel('Depth (m)', fontsize=12)
    plt.title(f'2D SH at {k*dt*1000:.2f} ms')
    k += mskip

anim = animation.FuncAnimation(fig, SH_2D, frames=int((nt - 2 * mskip) / mskip), interval=50)
anim.save('wave_propagation.mp4', writer='ffmpeg', fps=20)
