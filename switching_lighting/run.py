#Code for testing out how to use mpi4py
import subprocess #For interacting with the C++ code
import numpy as np

def simulation(parameter):
    params = ['%.16g' % a for a in parameter]
    exe_path = './Subprocess'
    subprocess.check_output([exe_path] + params)
N = 40
div = np.linspace(0.0,1.0,N)
out = np.zeros((N,3), dtype=float)

if div[len(div)-1] != -1.0:
    for i in range(len(div)):
	# Parameters for simulating ON
        p = [1.0*div[i]]
        simulation(p)
        F = np.loadtxt("SimulationDataForce.txt")
        C = np.loadtxt("SimulationDataCTot.txt")
        E = np.loadtxt("SimulationDataETot.txt")
        M = np.loadtxt("ModelOutputData.txt")
        if i==0:
            t = np.loadtxt("SimulationDataTime.txt")
            t_filename = f"SimulationDataTime_use_set1.txt"
            np.savetxt(t_filename, t, fmt='%16f')

        out[i,0] = min(F)
        out[i,1] = C[0]
        out[i,2] = E[0]

        # Format the filename with leading zeros
        F_filename = f"SimulationDataForce_{(i+1):02d}_set1.txt"
        C_filename = f"SimulationDataCTot_{(i+1):02d}_set1.txt"
        E_filename = f"SimulationDataETot_{(i+1):02d}_set1.txt"
        M_filename = f"ModelOutputData_{(i+1):02d}_set1.txt"
        out_file = "Output_E_C_F.txt"

        # Save the numpy array to a text file using np.savetxt
        np.savetxt(F_filename, F, fmt='%16f')
        np.savetxt(C_filename, C, fmt='%16f')
        np.savetxt(E_filename, E, fmt='%16f')
        np.savetxt(M_filename, M, fmt='%16f', delimiter='\t')
        np.savetxt(out_file, out, fmt='%16f', delimiter='\t')


else:
    p = [1.0]
    simulation(p)