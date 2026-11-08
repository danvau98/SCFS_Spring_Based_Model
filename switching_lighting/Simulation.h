#pragma once
//Carries and intialises the simulation
void simulation(double* t, double* SimDataL1, double* SimDataL2, double* SimDataL3, double* SimDataL4, double* SimDataBTot, double* SimDataCTot, double* SimDataETot, double* SimDataECTot, double* SimDataForce, double* SimDataTime, double* LT, double L1, double L2, double L3, double L4, double at, double bt, double dt, double TimeAPP, double pause, double div, int Nt, int NSim);
//Solves for the spring constants when all layers are connected
double* k_solve(double* k, double param_ECM, double param_UPPER, double param_cell, double param_LOWER, double B, double EC, double rho, double dx, int Nx, int tn);
//Solves for the spring constants when the lower integrin layer is not connected
double* k_solve_2(double* k, double param_ECM, double param_UPPER, double param_cell, double param_LOWER, double B, double EC, double rho, double dx, int Nx, int tn);
//Solves for the spring constants when the cell is no longer attached to 
//the AFM probe 
double* k_solve_3(double* k, double param_ECM, double param_UPPER, double param_cell, double param_LOWER, double B, double EC, double rho, double dx, int Nx, int tn);
//Generate the displacement for the AFM probe when the experimental data is inputted
double* L_T_solve(double* LT_d, double* t, int tn, int section, std::vector<std::vector<std::vector<double>>> ExpData);
//Generates the data for the AFM probe displacement when experimental data is included
double* LT_Set(double* LT, double* time, double L3, double L4, double dt, double speed, double pause, double timeApp, double timeRet, double initDist, int Nt);
//The rhs of the ODE for the integrin binding in the layer between the cell and AFM probe
double b_rhs(double d, double B, double eta, double** b, double t, int i, int tn);
//Solves the ODE for the inegrin binding between the cell and the AFM probe
double** solve_Diff_upper(double** b, double* x, double* t, double eta, double dx, double dt, double L, int Nx, int tn, int Nt);
//The right-hs of the ODE for the integrin binding to the biomaterial in the contracted state
double c_rhs(double d, double E_C, double gamma_bar, double theta_bar, double** c, double** e, int lambda1, int lambda2, double t, int i, int tn);
//The right-hs of the ODE for the integrin binding to the biomaterial in the extended state 
double e_rhs(double d, double E_C, double gamma_hat, double theta_hat, double** c, double** e, int lambda1, int lambda2, double t, int i, int tn);
//Solves the ODE for the integrin layer inbetween the cell and the biomaterial 
double** solve_Diff_lower(double** e, double** c, double* x, double* t, int lambda1, int lambda2, double dx, double dt, double L, double a, double b, double gamma_bar, double gamma_hat, double theta_bar, double theta_hat, int Nx, int tn, int Nt);
//Calculates the force in the same manner as the experimental data
double Force_calc(double k, double Ln, double Ln_d);