import numpy as np
from time import time
### TODO Figure out why ###
### (lattice[:,:] != lattice[i,:]).all(axis=1)] 
### ~(lattice[:,:] == lattice[i,:]).all(axis=1)) 


a = 5.26 # Angstrom
sigma =  3.4 #Angstrom
fcc_cell = a*np.array( [ [ 0, 0, 0 ], [ .5, .5, 0], [ 0.5, 0, .5], [ 0, 0.5, 0.5] ] )
kb = 8.617e-5
eps = 120 *kb#times kb

x = np.array (  [1, 0, 0,] )
y = np.array (  [0, 1, 0,] )
z = np.array (  [0, 0, 1,] )

def perform(func):

    def timer(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f' took {t2-t1} s')
        return result
    return timer


def translate_cell(a, x_i, y_i, z_i ):
	return fcc_cell + a*(x_i*x + y_i*y +z_i*z)

@perform
def genlattice( x_i, y_i, z_i, fhand):
    counter = 0
    for i in range(x_i):
        for j in range(y_i):
            for k in range(z_i):
                for el,n in enumerate(translate_cell(a, i, j, k)):
                    if el == 0: 
                        fhand.write('X\t')
                    else:
                        fhand.write('y\t')
                    for num in n:
                        fhand.write(f'{num}\t')
                    fhand.write(f'\n')
                    print(counter)
                    counter+=1

def lj_pot_vector( r_1, r_2, L,eps, sigma, cutoff= False):
    if cutoff:
        pass
    else:
        r = np.linalg.norm((r_2-r_1),axis=1)
        sigr = sigma/r
        sig6r = sigr*sigr*sigr*sigr*sigr*sigr
        sig12r = sig6r * sig6r
        return 4*eps*( sig12r - sig6r )

def lj_pot( r_1, r_2, L,eps, sigma, cutoff= True):
    if cutoff:
        r = (r_2-r_1)
        r = r - np.round( r/L ) * L
        r = np.linalg.norm(r)
        if r>L/2:
            return 0
        sigr = sigma/r
        sig6r = sigr*sigr*sigr*sigr*sigr*sigr
        sig12r = sig6r * sig6r
        return 4*eps*( sig12r - sig6r )
    else:
        r = np.linalg.norm((r_2-r_1))
        sigr = sigma/r
        sig6r = sigr*sigr*sigr*sigr*sigr*sigr
        sig12r = sig6r * sig6r
        return 4*eps*( sig12r - sig6r )

#@perform
def calcpotentials(lattice, sigma:float, eps:float,**kwargs):
    sitepots = np.zeros( lattice.shape[0] )
    for i in range(sitepots.shape[0] ):
        sitepots[i] += np.sum( lj_pot_vector( lattice[ i, : ], lattice[~(lattice[:,:] == lattice[i,:]).all(axis=1)], eps, sigma ) )
    return sitepots/2
    
#@perform
def easymode(lattice, L, sigma:float, eps:float,**kwargs):
    sitepots = np.zeros( lattice.shape[0] )
    for i in range(sitepots.shape[0] ):
        for j in range(i+1,sitepots.shape[0] ):
            if i==j:
                continue
            en = lj_pot(lattice[i], lattice[j], L,eps, sigma)
            sitepots[i]+= en
            sitepots[j]+= en
    return sitepots

if __name__ == '__main__':
    fname = "fcclattice256.xyz"

    with open(fname, 'r') as fhand:
        val =  np.sum( easymode( fhand, 4*a,sigma, eps ,test=True ), axis = 0 )
        print(f'The force per molecule: {val/(256)}' )

