import numpy as np
from relaxations import relax_hA
from parameters import *

# define objective function
def objective(params):
    '''
    Takes a given set of parameters and returns the total sum-squared error
    against target calues 
    '''

    # generate solution
    sol_jit = relax_hA(params) 

    # calculate error after equilibrium is reached
    T = np.arange(0,10000,0.1)
    T_eq = 9000
    i_eq = [i for i in range(len(T)) if T[i] >= T_eq]

    pow = np.array([s[6]*P for s in sol_jit])[i_eq]
    pow_error = sum((pow - P)**2)  

    fuel_in = np.array([s[15] for s in sol_jit])[i_eq]
    fuel_in_error = sum((fuel_in-F_to_K(1209))**2)

    fuel_out = np.array([s[1] for s in sol_jit])[i_eq]
    fuel_out_error = sum((fuel_out-F_to_K(1522))**2)

    coolant_in = np.array([s[25] for s in sol_jit])[i_eq]
    coolant_in_error = sum((coolant_in-F_to_K(1226))**2)

    coolant_out = np.array([s[4] for s in sol_jit])[i_eq]
    coolant_out_error = sum((coolant_out-F_to_K(1335))**2)

    tot = pow_error + fuel_in_error + fuel_out_error + coolant_in_error + coolant_out_error

    return tot
