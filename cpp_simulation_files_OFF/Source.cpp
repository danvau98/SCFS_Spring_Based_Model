//------------------------------------------------------------------------
//Rough Outline for how this script works
// - Setup and initialise the various variables and arrays
// - Import the experimental data that is going to be used in the 
// parameter estimation.
// - Setup LT_d on the experimental distance that the AFM tip has moved
// - Carry out an intial run of the simulation with the AFM distance data
// - When EC=0 the value for force is outputted and compared using a cost 
// function against the force that was found in the experiment
// - Optimise the cost function using the L-M algo so that 
//	
//------------------------------------------------------------------------
//This section will contain the libraries for the various algorithms that| 
//will be required to fit these parameters								 |
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "Reading_and_Writing.h"
#include "Useful.h"
#include "Simulation.h"
//------------------------------------------------------------------------
std::vector<std::vector<double>> parseFile(const std::string& filename);
int lambda1_val_loc1(int i);
int lambda2_val_loc1(int i);
int lambda1_val_loc2(int i);
int lambda2_val_loc2(int i);
int lambda1_val_loc3(int i);
int lambda2_val_loc3(int i);
int lambda1_val_loc4(int i);
int lambda2_val_loc4(int i);
void Sim_Selection(double* SimData, double* SimData_New, double* ExpData, int NExp, int NSim);
//------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    //This section is passing the experimental data into the script
    std::vector<std::vector<double>> ExpData;
    ExpData = parseFile("avg_OFF_all_locs.txt");
    std::string div_str(argv[5]);
    double div = std::stod(div_str);


    int i = 0;

    int lambda1;
    int lambda2;
    //OFF state
    if (i == 0)
    {
        lambda1 = 1; //OFF: 1, ON: 0 # change this
        lambda2 = 0; //OFF: 0, ON: 1 # change this
    }
    //ON state
    else
    {
        lambda1 = 0; //OFF: 1, ON: 0 # change this
        lambda2 = 1; //OFF: 0, ON: 1 # change this
    }
    //-------------------------------------------
    //Various variables and arrays that are
    //needed for parameter estimation 
    double pause = 1.0; //3.0;
    //pause = pause * div;
    //pause = pause;
    double speed = 30.0 * pow(10, -6.0) * pow(10.0, 9); //speed the AFM probe approaches and retracts
    speed = speed * div;
    double StartDistAPP;
    double TimeAPP, TimeRET, TimeTotal, dt, InitTime;
    int NExp, NExpRow, Nt, NSim; // NParam;

    NExp = ExpData.size();
    TimeAPP = ExpData[NExp - 1][0]; //Time taken for approach
    TimeRET = ExpData[NExp - 1][0]; //Time taken for retraction
    InitTime = 0.0;
    TimeTotal = TimeAPP + pause + TimeRET;
    dt = 1e-04;
    Nt = round(TimeTotal/dt);
    dt = (TimeTotal - InitTime) / ((double) Nt - 1.0);

    double* t = new double[Nt+1];
    double* LT = new double[Nt+1];

    double L1 = 100.0;
    double L2 = 10.0;
    double L3 = 200.0;
    double L4 = 10.0;

    //std::string param_str_ECM(argv[1]);
    //double param_ECM = std::stod(param_str_ECM);
    double param_ECM = 1e-01;

    std::string param_str_UPPER(argv[1]);
    double param_UPPER = std::stod(param_str_UPPER);
    //param_UPPER = param_UPPER * div;

    //std::string param_str_cell(argv[1]);
    //double param_cell = std::stod(param_str_cell);
    double param_cell = 5.0e-02;

    //std::string param_str_LOWER(argv[4]);
    //double param_LOWER = std::stod(param_str_LOWER);
    double param_LOWER = 1.0;

    //std::string param_str_alpha(argv[2]);
    //double param_alpha = std::stod(param_str_alpha);
    //param_alpha = param_alpha * div;
    double param_alpha = 1e-04;

    std::string param_str_beta(argv[2]);
    double param_beta = std::stod(param_str_beta);
    //param_beta = param_beta * div;

    //std::string param_str_gamma_bar(argv[4]);
    //double param_gamma_bar = std::stod(param_str_gamma_bar);
    double param_gamma_bar = 1e-04;

    std::string param_str_gamma_hat(argv[3]);
    double param_gamma_hat = std::stod(param_str_gamma_hat);

    //std::string param_str_theta_bar(argv[6]);
    //double param_theta_bar = std::stod(param_str_theta_bar);
    // param_theta_bar = param_theta_bar * div;
    double param_theta_bar = 1e-04;

    std::string param_str_theta_hat(argv[4]);
    double param_theta_hat = std::stod(param_str_theta_hat);
    // param_theta_hat = param_theta_hat * div;

    //This is for comparing the force data with experimental
    double* SimData_L1;
    double* SimData_L2;
    double* SimData_L3;
    double* SimData_L4;
    double* SimData_BTot;
    double* SimData_CTot;
    double* SimData_ETot;
    double* SimData_ECTot;
    double* SimData_Force;
    double* SimData_Time;
    double* SimData_New_L1;
    double* SimData_New_L2;
    double* SimData_New_L3;
    double* SimData_New_L4;
    double* SimData_New_BTot;
    double* SimData_New_CTot;
    double* SimData_New_ETot;
    double* SimData_New_ECTot;
    double* SimData_New_Force;
    double* SimData_New_Time;
    double* ExpDataRET;
    double* ExpDataRETDist;

    double tol = 1e-10;

    std::ofstream ExpFile;
    ExpFile.open("ExperimentalData.txt");
    std::ofstream ExpFileDist;
    ExpFileDist.open("ExperimentalDataDist.txt");
    std::ofstream SimFileL1;
    SimFileL1.open("SimulationDataL1.txt");
    std::ofstream SimFileL2;
    SimFileL2.open("SimulationDataL2.txt");
    std::ofstream SimFileL3;
    SimFileL3.open("SimulationDataL3.txt");
    std::ofstream SimFileL4;
    SimFileL4.open("SimulationDataL4.txt");
    std::ofstream SimFileBTot;
    SimFileBTot.open("SimulationDataBTot_OFF.txt");
    std::ofstream SimFileECTot;
    SimFileECTot.open("SimulationDataECTot_OFF.txt");
    std::ofstream SimFileCTot;
    SimFileCTot.open("SimulationDataCTot_OFF.txt");
    std::ofstream SimFileETot;
    SimFileETot.open("SimulationDataETot_OFF.txt");
    std::ofstream SimFileForce;
    SimFileForce.open("SimulationDataForce_OFF.txt");
    std::ofstream SimFileTime;
    SimFileTime.open("SimulationDataTime_OFF.txt");

    NExp = ExpData.size();

    TimeAPP = ExpData[NExp - 1][0]; //Time taken for approach
    TimeRET = ExpData[NExp - 1][0]; //Time taken for retraction
    TimeTotal = TimeAPP + pause + TimeRET;
    dt = (TimeTotal - InitTime) / ((double) Nt - 1.0);

    StartDistAPP = TimeAPP * speed; //This will be the final distnce based this is retraction data

    NSim = floor(TimeRET / dt);

    SimData_L1 = new double[NSim]; //For initial simulated data
    SimData_L2 = new double[NSim]; //For initial simulated data
    SimData_L3 = new double[NSim]; //For initial simulated data
    SimData_L4 = new double[NSim]; //For initial simulated data
    SimData_BTot = new double[NSim]; //For initial simulated data
    SimData_CTot = new double[NSim]; //For initial simulated data
    SimData_ETot = new double[NSim]; //For initial simulated data
    SimData_ECTot = new double[NSim]; //For initial simulated data
    SimData_Force = new double[NSim];
    SimData_Time = new double[NSim];

    SimData_New_L1 = new double[NExp]; //For simulated data after comparison to original 
    SimData_New_L2 = new double[NExp]; //For simulated data after comparison to original 
    SimData_New_L3 = new double[NExp]; //For simulated data after comparison to original 
    SimData_New_L4 = new double[NExp]; //For simulated data after comparison to original 
    SimData_New_BTot = new double[NExp]; //For simulated data after comparison to original
    SimData_New_CTot = new double[NExp]; //For simulated data after comparison to original
    SimData_New_ETot = new double[NExp]; //For simulated data after comparison to original
    SimData_New_ECTot = new double[NExp]; //For simulated data after comparison to original
    SimData_New_Force = new double[NExp]; //For simulated data after comparison to original
    SimData_New_Time = new double[NExp];

    ExpDataRET = new double[NExp];
    ExpDataRETDist = new double[NExp];

    //Assuming that there is zero acceleration in the AFM tip 
    //We now setup the time domain
    tDomain(t, InitTime, TimeTotal, dt, Nt);

    //Setting the LT based on the experimental data and time length
    LT_Set(LT, t, L3, L4, dt, speed, pause, TimeAPP, TimeRET, abs(StartDistAPP), Nt);

    simulation(t, SimData_L1, SimData_L2, SimData_L3, SimData_L4, SimData_BTot, SimData_CTot, SimData_ETot, SimData_ECTot, SimData_Force, SimData_Time, LT, L1, L2, L3, L4, InitTime, TimeTotal, dt, param_ECM, param_UPPER, param_cell, param_LOWER, param_alpha, param_beta, param_gamma_bar, param_gamma_hat, param_theta_bar, param_theta_hat, TimeAPP, pause, Nt, NSim, lambda1, lambda2);
                
    for (int j = 0; j < NExp; j++)
    {
        //Should pull out the force data and scale them based on the final force
        //value recorded
        ExpFile << (ExpData[j][1] * pow(10.0, 9.0)) << std::endl; 
        ExpDataRET[j] = (ExpData[j][1] * pow(10.0, 9.0)); 
    }
     
    //Now that both files are present the force values from the experiment need to
    //correspond with the simulations and this will be done in a separate function.
    Sim_Selection(SimData_L1, SimData_New_L1, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_L2, SimData_New_L2, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_L3, SimData_New_L3, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_L4, SimData_New_L4, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_BTot, SimData_New_BTot, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_CTot, SimData_New_CTot, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_ETot, SimData_New_ETot, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_ECTot, SimData_New_ECTot, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_Force, SimData_New_Force, ExpDataRET, NExp, NSim);
    Sim_Selection(SimData_Time, SimData_New_Time, ExpDataRET, NExp, NSim);
    //Now both the Experimental and coressponding data has to be exported to a textfile.
    double sum = 0;
    for (int j = 0; j < NExp; j++)
    {
        SimFileL1 << SimData_New_L1[j] << std::endl;
        SimFileL2 << SimData_New_L2[j] << std::endl;
        SimFileL3 << SimData_New_L3[j] << std::endl;
        SimFileL4 << SimData_New_L4[j] << std::endl;
        SimFileBTot << SimData_New_BTot[j] << std::endl;
        SimFileCTot << SimData_New_CTot[j] << std::endl;
        SimFileETot << SimData_New_ETot[j] << std::endl;
        SimFileECTot << SimData_New_ECTot[j] << std::endl;
        SimFileForce << SimData_New_Force[j] << std::endl;
        SimFileTime << SimData_New_Time[j] << std::endl;

        sum += pow((SimData_New_Force[j] - (ExpData[j][1] * pow(10.0, 9.0))),2.0);
    }
    std::cout << pow(sum,0.5) << std::endl;

    //Closing the textfile
    ExpFile.close();
    ExpFileDist.close();
    SimFileL1.close();
    SimFileL2.close();
    SimFileL3.close();
    SimFileL4.close();
    SimFileBTot.close();
    SimFileCTot.close();
    SimFileETot.close();
    SimFileECTot.close();
    SimFileForce.close();
    SimFileTime.close();

    //Deleting dynamically allocated memory
    delete[] t; 
    delete[] SimData_L1;
    delete[] SimData_L2;
    delete[] SimData_L3;
    delete[] SimData_L4;
    delete[] SimData_BTot;
    delete[] SimData_CTot;
    delete[] SimData_ETot;
    delete[] SimData_ECTot;
    delete[] SimData_Time;
    delete[] SimData_New_L1;
    delete[] SimData_New_L2;
    delete[] SimData_New_L3;
    delete[] SimData_New_L4;
    delete[] SimData_New_BTot;
    delete[] SimData_New_CTot;
    delete[] SimData_New_ETot;
    delete[] SimData_New_ECTot;
    delete[] SimData_New_Time;
    delete[] ExpDataRET;
    return 0;
}
//------------------------------------------------------------------------
std::vector<std::vector<double>> parseFile(const std::string& filename) 
{
    std::ifstream inputFile(filename);

    if (!inputFile) {
        std::cout << "Unable to open file";
        exit(1);
    }

    std::string line;
    std::vector<std::vector<double>> section;

    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::vector<double> row(2);
        for (int i = 0; i < 2; i++) 
        {
            iss >> row[i];
        }
        section.push_back(row);
    }

    inputFile.close();

    return section;
}
void Sim_Selection(double* SimData, double* SimData_New, double* ExpData, int NExp, int NSim)
{
    int n = ceil(NSim / NExp);
    int i = 0;
    for (int j = 0; j < NExp; j++)
    {
        SimData_New[j] = SimData[i];
        i = i + n;
    }
}
