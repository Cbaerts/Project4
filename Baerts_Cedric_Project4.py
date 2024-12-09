import numpy as np
import matplotlib.pyplot as plt


# method = input('Choose a numerical method, FTCS or Crank-Nicholson: ')
method = 'ftcs'
# number of grid points
nspace = 1000
# number of steps
ntime = 1000
# length of computation region
length = 200
# time increments
tau = 0.04

#Initiates position grid.
xArray = np.linspace(-length/2, length/2, nspace)

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
    b: float
        Value of one below the diagonal
    d: float
        Value of the diagonal
    a: float
        Value of one above the diagonal
    Return
    ---------
    matrix: ndarray
        The tridiagonal matrix based on the inputted values
    '''

    matrix = np.zeros((N, N))
    aMatrix, bMatrix, dMatrix = a*np.eye(N, k=1), b*np.eye(N, k=-1), d*np.eye(N, k=0)
    matrix = matrix + aMatrix + bMatrix + dMatrix

    matrix[0, -1] = b
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
    psi = np.exp(imaginary*k0*xArray)*np.exp(-(xArray-x0)**2/(2*sigma0**2))/(np.sqrt(sigma0*np.sqrt(np.pi)))    
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

def sch_eqn(nspace, ntime, tau, method='ftcs', length=200, wparam = [10, 0, 0.5]):
    '''
    nspace is the number of spatial grid points, 
    ntime is the number of time steps to be evolved,
    tau is the time step to be used
    method: string, either ftcs or crank
    length: float, size of spatial grid. Default to 200 (grid extends from -100 to +100)
    wparam: list of parameters for initial condition [sigma0, x0, k0]. Default [10, 0, 0.5].
    '''
    tArray = np.arange(0, ntime*tau, tau)
    wavePacket = make_initialcond(wparam, xArray)
    h = length/(nspace-1)
    psi = np.zeros((ntime, nspace), dtype=complex)
    psi[0, :] = wavePacket
    I = np.eye(nspace)
    constant = -hBar**2/(2*mass*h**2)

    if method.lower() in ['ftcs', '1']:
        # eqn 9.31
        H = make_tridiagonal(nspace, constant, 1-2*constant, constant)
        H = H * (imaginary*tau/hBar)
        eig = spectral_radius(I-H) 

        if eig > 1:
            print(eig)
            return 'This matrix is unstable'
        else: 
            print('This integration is gonna be cool')
            for n in range(ntime-1):
                psi[n+1, :] = np.dot((I-H), psi[n, :])

    elif method.lower() in ['crank', 'crank-nicholson', '2']:
        # 9.46
        H = make_tridiagonal(nspace, constant, -2*constant, constant)
        # 9.40
        HN = np.dot(np.linalg.inv(I + (imaginary*tau/2/hBar)*H), I - (imaginary*tau/2/hBar)*H)
        for n in range(1, ntime):
            psi[n, :] = np.dot(HN, psi[n-1, :])
            
    else:
        return 'Choose an actual method you donkey!'

    prob = np.abs(psi*np.conjugate(psi))
    return psi, xArray, tArray, prob

def sch_plot(solution, figure):
    psi, xArray, tArray, prob = solution
    if figure[0].upper() == 'Y':
        plt.subplot()
        for i,T in enumerate(tArray):
                if i%100 == 0:
                    plt.plot(xArray, np.real(psi[i]), label=f"Time = {T}")
    plt.title("Solution to The Schrodinger Equation")
    plt.legend()
    plt.xlabel("Position")
    plt.ylabel("Ψ")
    # plt.savefig("BaertsCedric_fig1_Project4")
    plt.show()
    plt.close()
    if figure[1].upper() == 'Y':
        for i, T in enumerate(tArray):
                if i%100 == 0:
                    plt.plot(xArray, np.real(prob[i]), label=f"Time = {T}")
    plt.title("Probabilty of The Schrodinger Equation")
    plt.legend()
    plt.xlabel("Position")
    plt.ylabel("Probablity %")
    # plt.savefig("BaertsCedric_fig2_Project4")
    plt.show()
    plt.close()
    
    return
# figure = [input("Do you want psi plot (Y/N): "), input("Do you want prob plot (Y/N): ")]
figure = ['Y','Y']
testy = sch_eqn(nspace, ntime, tau)
if isinstance(testy, str):
    print(testy)
else:
    apple = sch_plot(testy, figure)