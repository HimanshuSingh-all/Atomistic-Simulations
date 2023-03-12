import numpy as np
from fcclattice import easymode
from time import time
### TODO Figure out why ###
### (lattice[:,:] != lattice[i,:]).all(axis=1)] 
### ~(lattice[:,:] == lattice[i,:]).all(axis=1)) 

fname = "fcclattice256.xyz"
traj = 'trajdiff.txt'
macro = 'macrostate.txt'
WRITECNTR = 0
a = 5.26
sigma =  3.4
kb = 8.617e-5
eps = 120*kb #times kb
dt = 0.1
m = 39.948
_hdtbm = dt*0.5/m
_hdt2bm =_hdtbm * dt 
NSTEPS = 500


np.seterr(all = 'raise') 
def lj_force(r_1, r_2, L, eps, sigma):
    rvec = r_1 - r_2
    rvec = rvec - L * np.round_( rvec/L )
    r = np.linalg.norm(rvec)
    if r>L/2:
        return np.array([0, 0, 0])
    try:
        sigr = sigma/r
    except ZeroDivisionError: 
        print(r_1, r_2, rvec, r)
    except FloatingPointError: 
        print(r_1, r_2, rvec, r)

    sig6r = sigr*sigr*sigr*sigr*sigr*sigr
    sig12r = sig6r * sig6r
    return 4 * eps * ( - 12 * sig12r + 6 * sig6r ) * rvec/(r**2)

def lj_force_vec(r_1, r_2, L, eps, sigma):
    r =  (r_1-r_2)-L*np.round_((r_1-r_2)/L)
    normr = np.linalg.norm(r,axis=1)
    valid_indices = np.logical_and(normr<L/2, normr>0) 
    r = r[valid_indices]
    normr = normr[valid_indices]
    r2=normr**2
    sigr2 = sigma**2/r2
    sig6r = sigr2**3
    sig12r = sig6r * sig6r
    return (4* eps * ( - 12 * sig12r/(r2) + 6 * sig6r/(r2) ) * r.T).T


def calcforces(lattice, simboxlen, sigma:float, eps:float,**kwargs):
    L = simboxlen
    siteforces = np.zeros( (lattice.shape[0] , 3 ) )
    for i in range(siteforces.shape[0] ):
        siteforces[i,:] +=np.sum(lj_force_vec(lattice[i,:],lattice,4*a, eps, sigma) ) #np.sum( lj_force( lattice[ i, : ]*np.ones([lattice.shape[0]-1,lattice.shape[1]]), lattice[~(lattice[:,:] == lattice[i,:]).all(axis=1)], L ,eps, sigma ) )
    return siteforces
    
#@perform
def easyforce(lattiice, simboxlen, sigma:float, eps:float,**kwargs):
    L = simboxlen
    siteforces = np.zeros( (lattice.shape[0], 3) )
    for i in range(siteforces.shape[0]-1 ):
        for j in range(i+1,siteforces.shape[0] ):
            if i==j:
                continue
            try:
                en = lj_force(lattice[i], lattice[j], L, eps, sigma)
            except:
                print(i,j, lattice[i], lattice[j])
            siteforces[i]+= en
            siteforces[j]-= en
    #print("Easymode Done")
    print('sitforces: ',siteforces[0],siteforces[-1])
    return siteforces




def write_traj(traj, lattice, latticevels, forces):
    """
    Writes the trajectory to the file
    """
    if WRITECNTR:
        opmode = 'a+'
    else:
        opmode = 'w+' 

    with open( traj, opmode ) as traj:
        for i in range(lattice.shape[0]):
            for j in range(3):
                traj.write(f'{lattice[i,j]}\t')

            for j in range(3):
                traj.write(f'{latticevels[i,j]}\t')

            for j in range(3):
                traj.write(f'{forces[i,j]}\t')

            traj.write('\n')

def write_macrostate(macro, lattice, latticevels, Step):
    
    if WRITECNTR:
        opmode = 'a+'
    else:
        opmode = 'w+' 

    with open( macro, opmode ) as f:
        if not WRITECNTR:
            f.write(f'step\tk\tU\tT\n')
        K = 0.5 * m * np.sum( np.sum(latticevels**2) )
        U = np.sum( easymode(lattice, 4*a, sigma, eps) )
        T = 2/3* K /(lattice.shape[0]*kb)
        f.write(f'{Step}\t{K}\t{U}\t{T}\n')

    
def writevmd(fhand, traj):
    pass

lattice  = np.loadtxt( fname , skiprows = 2 , usecols = (1, 2 ,3 ) )
lattice[0] = np.array([1.2,1.1,0.0])
latticevels = np.zeros_like( lattice )
force = easyforce(lattice, 4*a, sigma, eps )
for i in range(NSTEPS):

    if i%1==0:
        write_traj(traj, lattice, latticevels, force)
        write_macrostate(macro, lattice, latticevels, i)
        WRITECNTR = 1
    print(i, latticevels[250])
    lattice = lattice + latticevels*dt +_hdt2bm*force
    # half step update of velocity
    latticevels = latticevels + force * _hdtbm
    #Calculate the force for the next step
    force = easyforce(lattice, 4*a, sigma, eps)
    latticevels = latticevels + force * _hdtbm

