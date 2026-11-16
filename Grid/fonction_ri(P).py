import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

#rayon de Bragg-Slater
rm = 

#le but est de plot notre fonction r_i
ANGSTROM_TO_BOHR = 1.8897259886
def r_i(P):
    i=np.arange(1,P+1)
    x_i=np.cos(i*np.pi/(P+1))
    r_i=rm*(1+x_i)/(1-x_i)
    return r_i

def plot_ri(P):
    r_is=r_i(P)
    i=np.arange(1,P+1)
    grid_bohr = r_is * ANGSTROM_TO_BOHR
    plt.figure(figsize=(8, 5))
    plt.plot(r_is, i, marker='o', markersize=4, linestyle='-', color='#1f77b4')
    plt.xlabel('Rayon r_i')
    plt.ylabel('indice i')
    plt.title(f"Distribution des poitns radiaux")
    plt.xscale('log')
    
    