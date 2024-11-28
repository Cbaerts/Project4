import numpy as np
import matplotlib.pyplot as plt


# method = input('Choose a numerical method, FTCS or Crank-Nicholson: ')
method = 'ftcs'
# nspace = int(input('Enter the number of grid points: '))
nspace = 200
# ntime = int(input('Enter number of steps: '))
ntime = 500
# length = float(input("Enter length of computation region: "))
length = 200
# tau = int(input("Enter a value for tau: "))
tau = 1

def make_initialcond(sigma0, k0,xArray):
    # From Lab 10. Liam "Handsome boy" Wilson wrote it. 
    """
    Takes a position array as input and returns the initial condition for the positions in the input array.
    The initial condition is a gaussian function times a cosine function.

    Args:
    sigma0 (float): Standard Deviation of the gaussian function.
    k0 (float): Frequency of cosine function.
    xArray (numpy array): Position grid

    Returns:
    a0 (Numpy Array): Initial condition corresponding to the input position array.
    """
    a0 = np.exp((-xArray**2)/(2*sigma0**2))*np.cos(k0*xArray)
    return a0

def spectral_radius(matrix):
    ''' 
    Calculates the eigenvalues of the inputted matrix

    Parameters
    ----------
    matrix: ndarray
        A NxN matrix that you want the absolute max eigenvalue for

    Return
    ---------
    np.max(abs(eigenvalues)): float
        The absolute max of the eigenvalues array
    '''
    eigenvalues = np.linalg.eig(matrix)[0]
    return np.max(abs(eigenvalues))

#Initiates position grid.
nSpace = 300
x_array = np.linspace(-length/2, length/2, nSpace)

#Parameters for initial condition function.
sigma = 0.2
k = 35
params = [sigma, k, x_array]


def sch_eqn(nspace, ntime, tau, method='ftcs', length=200, potential = [], wparam = [10, 0, 0.5]):
    sigma0, x0, k0 = wparam
    initCond = make_initialcond()
    if method == 'ftcs':
        # do stability check
        eig = spectral_radius()
        i = 1
        if eig >= tau:
            return 'This shit is unstable'
    elif method == 'crank':
        # Ignore stability check
        i = 1
    else:
        print('Choose an actual method you donkey!')
