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

    fig = plt.figure(figsize=(11.69, 8.27))

    # Reading standard stellar model
    m_sun, r_sun, T_sun, rho_sun, P_sun = np.loadtxt('bp2004stdmodel.dat', usecols=(0, 1, 2, 3, 4), skiprows = 23, unpack=True)
    rho_sun *= 1.e3
    P_sun *= 1.e-1
    
    # =======
    # Density
    # =======
    
    plt.subplot(2, 2, 1)
    plt.plot(r_sun, rho_sun, '--r', label = r'standard solar model')
    plt.plot(x, rho, 'b-', label = r'solar polytrope $(n = 3)$')
    plt.xlabel(r'$\frac{r}{R_\star}$')
    plt.ylabel(r'$\rho\,[\mathrm{kg}\,\mathrm{m}^{-3}]$')
    plt.yscale('log')
    plt.ylim(ymin = 1.e-1, ymax = 1.e6)
    plt.xlim(xmin = 0., xmax = 1.)
    plt.legend(fontsize = 'small')
    
    # ========
    # Pressure
    # ========

    plt.subplot(2, 2, 2)
    plt.plot(x, P, 'b-')
    plt.plot(r_sun, P_sun, '--r')
    plt.xlabel(r'$\frac{r}{R_\star}$')
    plt.ylabel(r'$P\,[\mathrm{N}\,\mathrm{m}^{-2}]$')
    plt.yscale('log')
    plt.ylim(ymin = 1.e5, ymax = 1.e17)
    plt.xlim(xmin = 0., xmax = 1.)
    
    # ===========
    # Temperature
    # ===========

    plt.subplot(2, 2, 3)
    plt.plot(x, T, 'b-')
    plt.plot(r_sun, T_sun, '--r')
    plt.xlabel(r'$\frac{r}{R_\star}$')
    plt.ylabel(r'$T\,[\mathrm{K}]$')
    plt.yscale('log')
    plt.ylim(ymin = 1.e4, ymax = 5.e7)
    plt.xlim(xmin = 0., xmax = 1.)

    # ====
    # Mass
    # ====

    plt.subplot(2, 2, 4)
    plt.plot(x, M/Msun, 'b-')
    plt.plot(r_sun, m_sun, '--r')
    plt.xlabel(r'$\frac{r}{R_\star}$')
    plt.ylabel(r'$\frac{m(r)}{M_{\odot}}$')
    #plt.yscale('log')
    plt.ylim(ymin = 0., ymax = 1.1)
    plt.xlim(xmin = 0., xmax = 1.)
    
    plt.savefig('pstar.pdf')
    plt.close()

