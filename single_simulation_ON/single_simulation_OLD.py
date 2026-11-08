# Code that runs a single simulation with a fixed set of parameter values
#-------------------------------------------------------------------------------
# Libraries that are used
import subprocess  # For interacting with the C++ code
import numpy as np  # For handling arrays
#-------------------------------------------------------------------------------
# Functions that are required
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

# Fixed parameter values (order: pu, beta, gamma2, theta2)
# Set these to the values you want to simulate with
p_fixed = np.array([0.186, 0.0150, 7.63, 1.54,0.015216272,4.782643566,1.436769815,3.000177,3.247279486])

# Run the single simulation
result = simulation(p_fixed)

# Getting the simulated force data and experimental data
fileForce = "SimulationDataForce_ON.txt"
fileExp = "ExperimentalData.txt"
Force_Exp = np.loadtxt(fileExp, dtype=float)
Force = np.loadtxt(fileForce, dtype=float)

# Computing the error against the experimental data
error = np.abs(Force - Force_Exp)
mean_error = error.mean()
print("Parameters:", p_fixed)
print("Objective value:", result)
print("Mean absolute error vs experiment:", mean_error)

# Exporting the parameter set and error into a text file
np.savetxt("Fixed_Param_Single_Run.txt", np.atleast_2d(p_fixed), delimiter='\t')
np.savetxt("Error_Single_Run.txt", np.array([mean_error]), delimiter='\t')