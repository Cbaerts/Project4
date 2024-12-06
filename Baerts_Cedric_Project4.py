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

#Initiates position grid.
nSpace = 300
xArray = np.linspace(-length/2, length/2, nSpace)

#Parameters for initial condition function.
sigma = 0.2
k = 35
params = [sigma, k, xArray]

# Constants
hBar = 1
mass = 1/2
imaginary = 1j

def make_tridiagonal(N, b, d, a):
    ''' 
    Generates a matrix as set by equation 9.48 in the text

    Parameters
    ----------
    N: Int
        The dimensions of the matrix NxN
    a: float
        Value of one above the diagonal
    b: float
        Value of one below the diagonal
    d: float
        Value of the diagonal

    Return
    ---------
    matrix: ndarray
        The tridiagonal matrix based on the inputted values
    '''

    matrix = np.zeros((N, N))
    aMatrix, bMatrix, dMatrix = a*np.eye(N, k=1), b*np.eye(N, k=-1), d*np.eye(N, k=0)
    matrix = matrix + aMatrix + bMatrix + dMatrix

    matrix[0, -1], matrix[-1, -2] = b, b
    matrix[-1, -1] = d
    matrix[-1, 0]  = a
    return matrix

def make_initialcond(wParam, xArray):
    """
    Takes a position array as input and returns the initial condition for the positions in the input array.
    The initial condition is a gaussian function times a cosine function.

    Args:
    sigma0 (float): Standard Deviation of the gaussian function.
    k0 (float): Average wave number
    xArray (numpy array): Position grid
    x0 (float): 

    Returns:
    a0 (Numpy Array): Initial condition corresponding to the input position array.
    """
    sigma0, x0, k0 = wParam
    psi = np.exp(imaginary*k0*xArray)*np.exp(-(xArray-x0)**2/(2*sigma0**2))/(sigma0*(np.pi)**0.5)**0.5
    return psi

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

def sch_eqn(nspace, ntime, tau, method='ftcs', length=200, potential = [], wparam = [10, 0, 0.5]):
    '''
    nspace is the number of spatial grid points, 
    ntime is the number of time steps to be evolved,
    tau is the time step to be used
    method: string, either ftcs or crank
    length: float, size of spatial grid. Default to 200 (grid extends from -100 to +100)
    potential: 1-D array giving the spatial index values at which the potential V(x) should be set to 1. Default to empty. For example, [25, 50] will set V[25] = V[50] = 1.
    wparam: list of parameters for initial condition [sigma0, x0, k0]. Default [10, 0, 0.5].
    '''

    wavePacket = make_initialcond(wparam, xArray)
    V = np.ones(nspace)
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
