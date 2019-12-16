#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', **{'family' : 'serif', 'serif' : ['Computer Modern Roman']})
plt.rc('text', usetex = True)

Msun = 1.9891e30 # kg
Rsun = 6.95508e8 # m

if __name__ == '__main__':
    data = 'pstar.txt'
    x, rho, P, T, M = np.loadtxt(data, usecols = (0, 1, 2, 3, 4), unpack = True)

    # fig = plt.figure(figsize=(11.69, 8.27))

    # Reading standard stellar model
    m_sun, r_sun, T_sun, rho_sun, P_sun = np.loadtxt('bp2004stdmodel.dat', usecols=(0, 1, 2, 3, 4), skiprows = 23, unpack=True)
    rho_sun *= 1.e3
    P_sun *= 1.e-1

    # xlabel = r'$m\,[M_{\odot}]$'
    xlabel = r'$m / M_\star$'
    
    # =======
    # Density
    # =======
    
    plt.plot(m_sun, rho_sun, '--r', label = r'standard solar model')
    plt.plot(M / Msun, rho, 'b-', label = r'solar polytrope $(n = 3)$')
    plt.xlabel(xlabel)
    plt.ylabel(r'$\rho\,[\mathrm{kg}\,\mathrm{m}^{-3}]$')
    plt.yscale('log')
    plt.ylim(ymin = 1.e-1, ymax = 1.e6)
    plt.xlim(xmin = 0., xmax = 1.)
    plt.legend()
    plt.savefig('pstar_density.pdf')
    plt.close()
    
    # ========
    # Pressure
    # ========

    plt.plot(M / Msun, P, 'b-', label = r'solar polytrope $(n = 3)$')
    plt.plot(m_sun, P_sun, '--r', label = r'standard solar model')
    plt.xlabel(xlabel)
    plt.ylabel(r'$P\,[\mathrm{N}\,\mathrm{m}^{-2}]$')
    plt.yscale('log')
    plt.ylim(ymin = 1.e11, ymax = 1.e17)
    plt.xlim(xmin = 0., xmax = 1.)
    plt.legend()
    plt.savefig('pstar_pressure.pdf')

    plt.close()
    
    # ===========
    # Temperature
    # ===========

    plt.plot(M / Msun, T, 'b-', label = r'solar polytrope $(n = 3)$')
    plt.plot(m_sun, T_sun, '--r', label = r'standard solar model')
    plt.xlabel(xlabel)
    plt.ylabel(r'$T\,[\mathrm{K}]$')
    plt.yscale('log')
    plt.ylim(ymin = 1.e5, ymax = 2.e7)
    plt.xlim(xmin = 0., xmax = 1.)
    plt.legend()
    plt.savefig('pstar_temperature.pdf')
    plt.close()
    

