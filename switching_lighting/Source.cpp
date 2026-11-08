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

    //-------------------------------------------
    std::string div_str(argv[1]);
    double div = std::stod(div_str);

    //Various variables and arrays that are
    //needed for parameter estimation
    double factor = 5.0-div;
 
    double pause = factor * 1.0; //3.0;
    double speed = 30.0 * pow(10, -6) * pow(10.0, 9); //speed the AFM probe approaches and retracts
    double StartDistAPP;
    double TimeAPP, TimeRET, TimeTotal, dt, InitTime;
    int NExp, NExpRow, Nt, NSim; // NParam;

    TimeAPP = 0.4; //Time taken for approach
    TimeRET = 0.4; //Time taken for retraction
    InitTime = 0.0;
    TimeTotal = TimeAPP + pause + TimeRET;
    dt = 5e-05;
    Nt = round(TimeTotal/dt);
    dt = (TimeTotal - InitTime) / ((double) Nt - 1.0);

    //NParam = 1;
    double* t = new double[Nt+1];
    double* LT = new double[Nt+1];

    double L1 = 100.0;
    double L2 = 10.0;
    double L3 = 200.0;
    double L4 = 10.0;

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
    SimFileBTot.open("SimulationDataBTot.txt");
    std::ofstream SimFileECTot;
    SimFileECTot.open("SimulationDataECTot.txt");
    std::ofstream SimFileCTot;
    SimFileCTot.open("SimulationDataCTot.txt");
    std::ofstream SimFileETot;
    SimFileETot.open("SimulationDataETot.txt");
    std::ofstream SimFileForce;
    SimFileForce.open("SimulationDataForce.txt");
    std::ofstream SimFileTime;
    SimFileTime.open("SimulationDataTime.txt");

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

    //Assuming that there is zero acceleration in the AFM tip 
    //We now setup the time domain
    tDomain(t, InitTime, TimeTotal, dt, Nt);

    //Setting the LT based on the experimental data and time length
    LT_Set(LT, t, L3, L4, dt, speed, pause, TimeAPP, TimeRET, abs(StartDistAPP), Nt);

    simulation(t, SimData_L1, SimData_L2, SimData_L3, SimData_L4, SimData_BTot, SimData_CTot, SimData_ETot, SimData_ECTot, SimData_Force, SimData_Time, LT, L1, L2, L3, L4, InitTime, TimeTotal, dt, TimeAPP, pause, div, Nt, NSim);

    //Now both the Experimental and coressponding data has to be exported to a textfile.
    double sum = 0;
    for (int j = 0; j < NSim; j++)
    {
        SimFileL1 << SimData_L1[j] << std::endl;
        SimFileL2 << SimData_L2[j] << std::endl;
        SimFileL3 << SimData_L3[j] << std::endl;
        SimFileL4 << SimData_L4[j] << std::endl;
        SimFileBTot << SimData_BTot[j] << std::endl;
        SimFileCTot << SimData_CTot[j] << std::endl;
        SimFileETot << SimData_ETot[j] << std::endl;
        SimFileECTot << SimData_ECTot[j] << std::endl;
        SimFileForce << SimData_Force[j] << std::endl;
        SimFileTime << SimData_Time[j] << std::endl;
    }

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
    return 0;
}
//------------------------------------------------------------------------
//Parssing the experimental data into an array that can be called on 
//and ouputted into arrays for estimation
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
//This function selects the appropriate points of simulated data based on
//the number of experimental data along the same time course.
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
