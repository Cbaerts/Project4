import numpy as np
import matplotlib.pyplot as plt


# method = input('Choose a numerical method, FTCS or Crank-Nicholson: ')
method = 'ftcs'
# nspace = int(input('Enter the number of grid points: '))
nspace = 300
# ntime = int(input('Enter number of steps: '))
ntime = 500
# length = float(input("Enter length of computation region: "))
length = 200
# C = float(input("Enter wave speed: "))
C = 1
# tau = int(input("Enter a value for tau: "))
tau = 1

def sch_eqn(nspace, ntime, tau, method='ftcs', length=200, potential = [], wparam = [10, 0, 0.5]):
    if method == 'ftcs':
        # do stability check
        i = 1
    elif method == 'crank':
        # Ignore stability check
        i = 1
    else:
        print('Choose an actual method you donkey!')