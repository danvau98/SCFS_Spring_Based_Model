# Code that should carry out the parameter estimation for the parameters pu,
# alpha, beta, gamma1, gamma2, theta1, theta2. Using the least squares method
#-------------------------------------------------------------------------------
# Libraries that are used
from scipy.stats import qmc
import subprocess #For interacting with the C++ code
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import least_squares
import numpy as np #For handling arrays
#-------------------------------------------------------------------------------
# Functions that are required
def my_cost_function(parameters):
    #Carrying out the simulation
    simulation(parameters)
    print(parameters)
    # Getting the simulated force data and experimental data
    fileForce = "SimulationDataForce.txt"
    fileExp = "ExperimentalData.txt"
    Force_Exp = np.loadtxt(fileExp, dtype=float)
    Force = np.loadtxt(fileForce, dtype=float)
    print(np.abs(Force-Force_Exp).mean())
    return np.abs(Force-Force_Exp)

def simulation(parameter):
    params = ['%.16g' % a for a in parameter]

    div = 1.0  # Replace with your intended scaling factor

    exe_path = './Subprocess'

    output = subprocess.check_output(
        [exe_path] + params + [str(div)]
    )
    print(params)
    print(float(output))
    return float(output)
#-------------------------------------------------------------------------------
# Main body of code
# Generating the parameter space using the LHS sampling
Num_params = 50
params = 4
error_cols = 1
# Arrays for outputting
est_p = np.zeros((Num_params,params))
rand_init_p = np.zeros((Num_params,params))
error_output = np.zeros((Num_params,error_cols))

# Main section of code lives here
sampler = qmc.LatinHypercube(d=4)
sample = sampler.random(n=Num_params)
# the lower and upper bound of the parameter space
#           pu | beta | g2 | h2
l_bounds = [0.001, 0.001, 0.01, 0.01]
u_bounds = [0.1, 0.1, 2.0, 2.0]
sample_scaled = qmc.scale(sample, l_bounds, u_bounds)

b = [0, 100]
bnds = Bounds(0.0, 100.0)

full_param_list = 4
for j in range(len(sample_scaled)):
    p0 = np.zeros(full_param_list)
    for i in range(len(p0)):
        p0[i] = sample_scaled[j,i]
        rand_init_p[j][i] = p0[i] 
        
    # Performing the estimation 
    p_new = minimize(simulation, p0, bounds = bnds, method = 'Nelder-Mead') #, method='lm', bounds = bnds

    # Outputting the estimated parameters
    for i in range(len(p_new.x)):
        est_p[j][i] = p_new.x[i]

    # Outputting the error from the estimation
    error_output[j] = (p_new.fun)

    print(j)

# Exporting the parameter and error array into a textfile
np.savetxt("Error_all_OFF_LHS_Nelder_L2_50.txt", error_output, delimiter='\t') # change this from location to location and state to state
np.savetxt("Est_Param_all_OFF_LHS_Nelder_L2_50.txt", est_p, delimiter='\t') # change this from location to location and state to state
np.savetxt("init_Param_all_OFF_LHS_Nelder_L2_50.txt", rand_init_p, delimiter='\t') # change this from location to location and state to state

