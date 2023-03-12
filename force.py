import numpy as np
import matplotlib.pyplot as plt
from time import time

dt, m, a = 0.5, 39.948, 5.26
kb = 8.617e-5
eps, sigma = 120*kb, 3.4
_hdbtm = 0.5*dt/m
_hdbt2m = _hdbtm*dt
NSTEPS = 5000
hm = 0.5*m
WP = 50 

init = 'fcclattice256.xyz'
lattice = np.loadtxt( init , skiprows = 2, usecols = (1,2,3) )
latticevels = np.zeros_like(lattice)


KE_to_Temp = 2 / (lattice.shape[0] * kb * 3)

def ljforce(r_1, r_2, boxlen, cutoff,**kwargs):
    
    r = r_1 - r_2
    r = r - boxlen*np.round(r/boxlen)
    normr = np.linalg.norm(r) 
    if np.linalg.norm(r) >= cutoff:
        return np.array( [0.0, 0.0, 0.0] )
    else:
        try:
            sigr = sigma/normr
            sig2r = sigr**2
            sig6r = sig2r**3
            sig12r = sig6r**2
            return  4 * eps * ( 12 * sig12r - 6 * sig6r ) * r/normr**2
        except FloatingPointError:
            print("before ",r_1, r_2, r_1 - r_2, np.round((r_1-r_2)/boxlen), boxlen * np.round((r_1-r_2)/boxlen) )
            print("after ",r, normr)
            print(kwargs['i'],kwargs['j'],':::')

np.seterr(all='raise')
def ljpot(r_1, r_2, boxlen, cutoff):
    
    r = r_1 - r_2
    r = r - boxlen*np.round(r/boxlen) 
    normr = np.linalg.norm(r) 
    if np.linalg.norm(r) >= cutoff:
        return 0.0
    else:
        try:
            sigr = sigma/normr
            sig2r = sigr**2
            sig6r = sig2r**3
            sig12r = sig6r**2
        except FloatingPointError:
            print("before ",r_1, r_2, r_1 - r_2)
            print("after ",r, normr)
        return  4 * eps * (  sig12r -  sig6r )

def calcforces(lattice):

    forces = np.zeros_like( lattice )
    for i, r_1 in enumerate(lattice[:]):
        for j in range(i+1, lattice.shape[0]):
            r_2 = lattice[j]
            force = ljforce(r_1, r_2, 4*a, 2*a, i=i, j=j)
            #print(f' pair ({i,j})')
            forces[i] += force
            forces[j] -= force
            
    return forces

def calcpot(lattice):

    pots = np.zeros( [lattice.shape[0]] )
    for i, r_1 in enumerate(lattice[:]):
        for j in range(i+1, lattice.shape[0]):
            r_2 = lattice[j]
            pot = ljpot(r_1, r_2, 4*a, 2*a)
            #print(f' pair ({i,j})')
            pots[i] += pot
            pots[j] += pot
    return pots
ti = time()            
lattice[0] = np.array([1.2 , 1.1, 0])
forces = calcforces(lattice)
with open('macroschange.txt','w+') as f:
    f.write(f' step\tK\tU\tT\n')
    for i in range(NSTEPS):
        if i%WP ==0:
            K = hm * np.sum( np.sum( latticevels**2 ) ) 
            U = np.sum(calcpot( lattice ))
            T = KE_to_Temp* K
            f.write(f'{i}\t{K:.6f}\t{U:.15f}\t{T:.6f}\n')
            #print(f'{i}\t{K:.6f}\t{U:.15f}\t{T:.6f}\n')

        latticevels = latticevels + _hdbtm*forces
        lattice = lattice + latticevels * dt
        forces = calcforces(lattice)
        latticevels = latticevels + _hdbtm * forces
        
    plt.legend()
    plt.show()
    tf = time()
print(f'TIme Taken: {tf-ti}')
