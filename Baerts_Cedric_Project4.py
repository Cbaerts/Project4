import numpy as np
import matplotlib.pyplot as plt

# https://github.com/Cbaerts/Project4

# Constants
hBar = 1
mass = 1/2
i = 1j

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
    matrix[-1, 0] = a
    return matrix

def make_initialcond(wParam, xArray):
    """
    Takes a position array as input and returns the initial condition for the positions in the input array.

    Parameters
    ----------
    wParam: nd.array
        [sigma0, x0, k0], where
        sigma0: Standard Deviation of the gaussian function.
        k0: Average wave number
        x0: Central position 
    xArray: nd.array 
        Position grid

    Returns
    ---------
    a0: nd.array 
        Initial condition for psi for the schrodinger equation.
    """
    sigma0, x0, k0 = wParam
    psi = np.exp(i*k0*xArray)*np.exp(-(xArray-x0)**2/(2*sigma0**2))/(np.sqrt(sigma0*np.sqrt(np.pi)))    
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
    eigenvalues: float
        The absolute max of the eigenvalues array
    '''
    eigenvalues = np.max(abs(np.linalg.eig(matrix)[0]))
    return eigenvalues

def sch_eqn(nspace, ntime, tau, method, length=200, wparam=[10, 0, 0.5]):
    '''
    Parameters
    ---------
    nspacet: int
        The number of spatial grid points, 
    ntime: int
        The number of time steps to be evolved,
    tau: int
        The time step to be used
    method: str
        Either ftcs or crank-nicholson
    length: float 
        size of spatial grid. Default to 200 (grid extends from -100 to +100)
    wparam: nd.array
        List of parameters for initial condition [sigma0, x0, k0]. Default [10, 0, 0.5].

    Return
    ---------
    psi: nd.array
        [ntime x nspace] shaped array which show the time evolution of psi
    xArray: nd.array
        1-D Array of all the points taken for psi
    tArray: nd.array
        1-D Array for all the time steps
    prob: nd.array
        [ntime x nspace] Array for the probabilities at all points and timesteps
    '''
    # Sets time stamps
    tArray = np.arange(0, ntime*tau, tau)
    # Create intial wavepacket
    wavePacket = make_initialcond(wparam, xArray)
    # Create empty complex array 
    psi = np.zeros((ntime, nspace), dtype=complex)
    # Set first line to the initial wavepacket
    psi[0, :] = wavePacket
    # Step distance
    h = length/(nspace-1)
    # Identity matrix
    I = np.eye(nspace)
    # Constant multiplied by both mehtods
    constant = -hBar**2/(2*mass*h**2)
    # Checks Method for FTCS
    if method.lower() in ['ftcs', '1']:
        # Makes Eqn 9.31 in Textbook
        H = make_tridiagonal(nspace, constant, 1-2*constant, constant)
        H = H * (i*tau/hBar)
        # Checks if the solution will be stable
        eig = spectral_radius(I-H) 
        # If solution is found to be unstable it will stop integration.
        if eig > 1:
            return 'The solution will be unstable. \n I WILL NOT INTEGRATE'
        # If solution is stable intergation will begin
        else: 
            print('The solution will be stable')
            for n in range(ntime-1):
                # Equation 9.32 in textbook
                psi[n+1, :] = np.dot((I-H), psi[n, :])
    # Checks Method for CRANK
    elif method.lower() in ['crank', 'crank-nicholson', '2']:
        # Equation 9.46 in the Textbook
        H = make_tridiagonal(nspace, constant, -2*constant, constant)
        # Equation 9.40 in the Textbook
        HN = np.dot(np.linalg.inv(I + (i*tau/2/hBar)*H), I - (i*tau/2/hBar)*H)
        for n in range(1, ntime):
            # Iterates through 
            psi[n, :] = np.dot(HN, psi[n-1, :])    
    # Returns a statement if no actual method was selected  
    else:
        return 'Choose an actual method!'
    # Calculates probability
    # Under Equation 9.42 in the Textbook
    prob = np.abs(psi*np.conjugate(psi))
    return psi, xArray, tArray, prob

def sch_plot(solution, figure):
    '''
    Parameters
    ---------
    Solution: tuple
        All the returns from the sch_eqn function. 
        [psi, xArray, tArray, prob]
    figure: list
        2 Item list which states if you want to get a psi plot and probablility plot
    Return
    ---------
    '''
    # Set solutions to their respective arrays
    psi, xArray, tArray, prob = solution
    # Set a number so always 5 timesteps are shown
    num = int(len(psi)/5)
    # Checks if psi plot wants to be created
    # Plots psi
    if figure[0].upper() == 'Y':
        plt.subplot()
        for i,T in enumerate(tArray):
                if i % num == 0:
                    plt.plot(xArray, np.real(psi[i]), label=f"Time step= {i}")
        plt.plot(xArray, np.real(psi[-1]), label=f"Time step = {len(tArray)}")
    plt.title("Wave Function of the Schrodinger Equation")
    plt.legend()
    plt.xlabel("Position")
    plt.ylabel("Ψ (x)")
    # plt.savefig("BaertsCedric_fig1_Project4")
    plt.show()
    plt.close()
    # Checks if probability plot wants to be created
    # Plots probability
    if figure[1].upper() == 'Y':
        for i, T in enumerate(tArray):
                if i % num == 0:
                    plt.plot(xArray, np.real(prob[i]), label=f"Time step = {i}")
        plt.plot(xArray, np.real(prob[-1]), label=f"Time step = {len(tArray)}")
    plt.title("Probabilty of The Schrodinger Equation")
    plt.legend()
    plt.xlabel("Position")
    plt.ylabel("Probablity (%)")
    # plt.savefig("BaertsCedric_fig2_Project4")
    plt.show()
    plt.close()
    return

# Ask for input parameters
method = input('Choose a numerical method, 1) FTCS or 2) Crank-Nicholson: ')
nspace = int(input('Enter the number of grid points: '))
ntime = int(input('Enter number of steps: '))
length = float(input("Enter length of computation region: "))
tau = float(input("Enter a value for tau: "))
# Working FTCS parameters
# nspace = 500
# ntime = 500
# length = 200
# tau = 0.01

# Initiates position grid.
xArray = np.linspace(-length / 2, length / 2, nspace)
# Initial wave packet parameters
parameters = [10, 0, 0.5]

# Call your function with the provided parameters
solution = sch_eqn(nspace, ntime, tau, method, length, parameters)

# Checks if your solution returned an error or was right
if isinstance(solution, str):
    print(solution)
else:
# If an actual solution asks if you want to see the Psi and/or probability plot
    figure = [input("Do you want psi plot (Y/N): "), input("Do you want probability plot (Y/N): ")]
    plot = sch_plot(solution, figure)

