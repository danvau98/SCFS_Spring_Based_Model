#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <cstring>
#include <cctype>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include "Deformation.h"
#include "Rate_param.h"
#include "Useful.h"
#include "Simulation.h"
//Carries and intialises the simulation
void simulation(double* t, double* SimDataL1, double* SimDataL2, double* SimDataL3, double* SimDataL4, double* SimDataBTot, double* SimDataCTot, double* SimDataETot, double* SimDataECTot, double* SimDataForce, double* SimDataTime, double* LT, double L1, double L2, double L3, double L4, double at, double bt, double dt, double param_ECM, double param_UPPER, double param_cell, double param_LOWER, double alpha, double beta, double gamma_bar, double gamma_hat, double theta_bar, double theta_hat, double l1, double l2, double TimeAPP, double pause, int Nt, int NSim, int lambda1, int lambda2)
{
	std::ofstream ModelOutput;
	ModelOutput.open("ModelOutputData.txt");

	std::ofstream c_output_init;
	c_output_init.open("c_output_init.txt");

	std::ofstream e_output_init;
	e_output_init.open("e_output_init.txt");

	std::ofstream b_output_init;
	b_output_init.open("b_output_init.txt");

	double rho;
	double eta;
	rho = 5.0; 
	eta = 0.6;

	//----------------------------------------------|
	// Number of time points - Nt					|
	// Number of spatial points - Nx				|
	// Number of displaced layers - N				|
	// Initialising the initial and final condtions |
	// initial time - at							|
	// final time - bt								|
	// initial space - ax							|
	// final space - bx								|
	// initial ECM displacement - L1				|
	// initial toplayer integrin displacement - L2	|
	// initial cell displacement - L3				|
	// inital underlayer integrin displacement - L4	|
	// initial AFM displacement - LT				|
	// difference between integrin displacements - L|
	// The total number of bound integrins - BTot	|
	//----------------------------------------------|
	int Nx = 400;
	int N = 4;
	double ax = -250.0;
	double bx = 250.0;
	double LT_eq = L1 + L2 + L3 + L4;
	double L1_d, L2_d, L3_d, L4_d, LT_d, L, BTot, ECTot, ECTot_n_scale, CTot, ETot;
	double Force;
	L1_d = L1;
	L2_d = L2;
	L3_d = L3;
	L4_d = L4;
	LT_d = LT[0];

	//----------------------------------------------|
	// Initialising the domains						|
	// temporal - t									|
	// spatial - x									|
	// No. time points - Nt							|
	// No. spacial points - Nx						|
	// time step - dt								|
	// space step - dx								|
	//----------------------------------------------|
	double dx;
	double* x = new double[Nx];
	double* k = new double[N];
	double* k2 = new double[3];
	double* k3 = new double[3];
	double* k4 = new double[2];
	double* B_TOT = new double[Nt + 1];
	double* EC_TOT = new double[Nt + 1];
	dx = (bx - ax) / ((double)Nx - 1.0);
	xDomain(x, ax, bx, dx, Nx);

	//----------------------------------------------|
	// The 2D array that will be populated with the |
	// number of bound integrins at a given point x |
	// and a givent time t.							|
	// bound integrins toplayer - b(x,t)			|
	// bound integrins extended - e(x,t)			|
	// bound integrins contracted - c(x,t)			|												
	//----------------------------------------------|
	double** b = new double* [Nx];
	double** e = new double* [Nx];
	double** c = new double* [Nx];


	//----------------------------------------------|
	// This is the half wave period to create a  	|
	// large enough displacement to get the			| 
	// neccessary change.							|
	//----------------------------------------------|
	k4[0] = param_ECM;
	k4[1] = param_cell;

	//----------------------------------------------|
	// The main while loop that performs the		|
	// simulation over the desired time course		|
	// following the plan listed in the heading of	|
	// this file.									|
	//----------------------------------------------|
	double* vec1 = new double[Nx];
	double* vec2 = new double[Nx];
	double* vec3 = new double[Nx];
	double* vec4 = new double[Nx];
	double* vec5 = new double[Nx];
	int tn, l;
	for (int i = 0; i < Nt; i++)
	{
		B_TOT[i] = 0.0;
		EC_TOT[i] = 0.0;
	}

	for (int i = 0; i < Nx; i++)
	{
		b[i] = new double[Nt];
		e[i] = new double[Nt];
		c[i] = new double[Nt];
	}
	//Setting up the initial conditions for the binding 
	for (int i = 0; i < Nx; i++)
	{
		b[i][0] = 0.08/100.0;
		e[i][0] = 0.0;
		c[i][0] = 0.0;
	}
	tn = 0;
	l = 0;
	L4_d = LT_d - L1 - L2 - L3 - L4;
	while (t[tn + 1] <= bt)
	{
		//----------------------------------------------|
		// Following the steps outlined in the header of|
		// this file to solve the necessary equations in| 
		// the model.									|
		// Folowing the steps							|
		//  1. Caluclate the spring constants			|
		//	2. Solve Micro-ODEs							|
		//	3. Calculate the springs deforming			|
		//	4. Repeat									|
		//----------------------------------------------|

		//----------------------------------------------|
		//Calculating BTot for the total number of bound|
		// integrins across the spatial domain			| 
		//----------------------------------------------|
		for (int j = 0; j < Nx; j++)
		{
			vec1[j] = b[j][tn];
			vec2[j] = (beta * e[j][tn]) + (alpha * c[j][tn]) ;
			vec3[j] = (c[j][tn]) + (e[j][tn]);
			vec4[j] = c[j][tn];
			vec5[j] = e[j][tn];
		}
		BTot = trapz(vec1, dx, Nx);
		ECTot = trapz(vec2, dx, Nx);
		ECTot_n_scale = trapz(vec3, dx, Nx);
		CTot = trapz(vec4, dx, Nx);
		ETot = trapz(vec5, dx, Nx);

		//Calculating the latest spring constants.
		k_solve(k, param_ECM, param_UPPER, param_cell, param_LOWER, BTot, ECTot, rho, dx, Nx, tn);
		k_solve_2(k2, param_ECM, param_UPPER, param_cell, param_LOWER, BTot, ECTot, rho, dx, Nx, tn);
		k_solve_3(k3, param_ECM, param_UPPER, param_cell, param_LOWER, BTot, ECTot, rho, dx, Nx, tn);

		//Calculating the new deformation in LT.
		LT_d = LT[tn];

		Force = Force_calc(k[0], L1, L1_d);

		ModelOutput << t[tn] << "\t" << BTot << "\t" << ECTot << "\t" << Force << "\t" << LT_d << "\t" << L1_d << "\t" << L2_d << "\t" << L3_d << "\t" << L4_d << std::endl;

		//Calculating the various displacements
		if (BTot > 0.0) //Binding with the tip
		{
			if (ECTot > 1e-10) //Binding with the biomaterial
			{
				if (LT_d > LT_eq)
				{
					L1_d = L_master(t, k, LT_eq, LT_d, L1, 0, N, tn);
					L2_d = L_master(t, k, LT_eq, LT_d, L2, 1, N, tn);
					L3_d = L_master(t, k, LT_eq, LT_d, L3, 2, N, tn);
					L4_d = L_master(t, k, LT_eq, LT_d, L4, 3, N, tn);
				}
				else
				{
					L2_d = L2;
					L4_d = L4;
					L3_d = L_master(t, k, LT_eq, LT_d, L3, 2, N, tn);
					L1_d = L_master(t, k, LT_eq, LT_d, L1, 0, N, tn);
				}
			}
			else //No binding with the biomaterial 
			{
				if ((LT_d) > L4_d)
				{
					L1_d = L1;
					L2_d = L2;
					L3_d = L3;
					L4_d = LT_d - L1_d - L2_d - L3_d;
				}
				else //Three spring model
				{
					L2_d = L2;
					L3_d = L_master(t, k2, LT_eq, LT_d, L3, 1, 3, tn);
					L1_d = L_master(t, k2, LT_eq, LT_d, L1, 0, 3, tn);
					L4_d = LT_d - L2_d - L3_d - L1_d;
				}
			}
		}
		else //No binding with the tip
		{
			if (LT_d - L1_d - L2_d > L3_d) //No compression
			{
				if (ECTot_n_scale > 1e-09) //Binding in the lower layer
				{
					//Three spring model, neither the upper or the
					//lower are rigid
					L1_d = L_master(t, k3, LT_eq, LT_d, L1, 0, 3, tn);
					L3_d = L_master(t, k3, LT_eq, LT_d, L3, 1, 3, tn);
					L4_d = L_master(t, k3, LT_eq, LT_d, L4, 2, 3, tn);
					L2_d = LT_d - L1_d - L3_d - L4_d;
				}
				else //No binding in the lower layer or in the upper layer
				{
					//No binding in both the upper and lower layer
					//Thus giving the two spring model
					L3_d = L_master(t, k4, LT_eq, LT_d, L3, 1, 2, tn);
					L1_d = L_master(t, k4, LT_eq, LT_d, L1, 0, 2, tn);
					L2_d = LT_d - L1_d - L3_d;
					L4_d = LT_d - L1_d - L3_d;
				}
			}
			else //compression
			{
				if (ECTot_n_scale > 1e-09) //Bindidng in the lower layer
				{
					//springs with rigid integrins
					L4_d = L4;
					L1_d = L_master(t, k3, LT_eq, LT_d, L1, 0, 3, tn);
					L3_d = L_master(t, k3, LT_eq, LT_d, L3, 1, 3, tn);
					L2_d = LT_d - L4_d - L1_d - L3_d;
				}
				else //No binding in the lower layer
				{
					//two spring model with no binding in both
					//the upper and lower layer
					L3_d = L_master(t, k4, LT_eq, LT_d, L3, 1, 2, tn);
					L1_d = L_master(t, k4, LT_eq, LT_d, L1, 0, 2, tn);
					L2_d = LT_d - L1_d - L3_d;
					L4_d = LT_d - L1_d - L3_d;
				}
			}
		}
		L = L2_d - L2;
		solve_Diff_upper(b, x, t, eta, dx, dt, L, Nx, tn, Nt);
		L = L4_d - L4;
		solve_Diff_lower(e, c, x, t, lambda1, lambda2, dx, dt, L, alpha, beta, gamma_bar, gamma_hat, theta_bar, theta_hat, l1, l2, Nx, tn, Nt);
		B_TOT[tn] = BTot;
		EC_TOT[tn] = ECTot;
		//--------------------------------------------
		//Saving the simulated data to the array
		if (t[tn] >= (TimeAPP + pause))
		{
			// Write the whole spatial profile of c as one row (tab-separated)
    		for (int i = 0; i < Nx; i++)
    		{
        		c_output_init << c[i][tn];
        		if (i < Nx - 1) c_output_init << "\t";  // separator between columns
    		}
    		c_output_init << std::endl;  // new line = new time point

    		// Same for e if you still want it
    		for (int i = 0; i < Nx; i++)
    		{
        		e_output_init << e[i][tn];
        		if (i < Nx - 1) e_output_init << "\t";
    		}
    		e_output_init << std::endl;

			// Same for b if you still want it
    		for (int i = 0; i < Nx; i++)
    		{
        		b_output_init << b[i][tn];
        		if (i < Nx - 1) b_output_init << "\t";
    		}
    		b_output_init << std::endl;

			SimDataL1[l] = L1_d;
			SimDataL2[l] = L2_d;
			SimDataL3[l] = L3_d;
			SimDataL4[l] = L4_d;
			SimDataBTot[l] = BTot;
			SimDataCTot[l] = CTot;
			SimDataETot[l] = ETot;
			SimDataECTot[l] = ECTot_n_scale;
			SimDataForce[l] = Force;
			SimDataTime[l] = t[tn];
			l += 1;
		}
		//--------------------------------------------
			
		//Adjusting the new time one step forward 
		tn += 1;
	}
	//----------------------------------------------|
	// Deleting dynamically allocated memory		|
	//----------------------------------------------|
	destroyTwoDimenArray_double(b, Nx, Nt);
	destroyTwoDimenArray_double(e, Nx, Nt);
	destroyTwoDimenArray_double(c, Nx, Nt);

	delete[] vec1;
	delete[] vec2;
	delete[] vec3;
	delete[] vec4;
	delete[] B_TOT;
	delete[] EC_TOT;
	delete[] k;
	delete[] k2;
	delete[] k3;
	delete[] k4;

	ModelOutput.close();
	c_output_init.close();
	e_output_init.close();
	b_output_init.close();
}
//Solves for the spring constants when all layers are connected
double* k_solve(double* k, double param_ECM, double param_UPPER, double param_cell, double param_LOWER, double B, double EC, double rho, double dx, int Nx, int tn)
{
	/*This function return an array that contains the necessary spring
	constants for the various different displacements.*/
	k[0] = param_ECM;
	k[1] = param_UPPER * B; 
	k[2] = param_cell;
	k[3] = param_LOWER * (EC); 

	return k;
}
//Solves for the spring constants when the lower integrin layer is not connected
double* k_solve_2(double* k, double param_ECM, double param_UPPER, double param_cell, double param_LOWER, double B, double EC, double rho, double dx, int Nx, int tn)
{
	/*This function return an array that contains the necessary spring
	constants for the various different displacements.*/
	k[0] = param_ECM; 
	k[1] = param_UPPER * B;
	k[2] = param_cell;

	return k;
}
//Solves for the spring constants when the cell is no longer attached to 
//the AFM probe 
double* k_solve_3(double* k, double param_ECM, double param_UPPER, double param_cell, double param_LOWER, double B, double EC, double rho, double dx, int Nx, int tn)
{
	/*This function return an array that contains the necessary spring
	constants for the various different displacements.*/
	k[0] = param_ECM; 
	k[1] = param_cell;
	k[2] = param_LOWER * (EC); 

	return k;
}
//Generate the displacement for the AFM probe when the experimental data is inputted
double* L_T_solve(double* LT_d, double* t, int tn, int section, std::vector<std::vector<std::vector<double>>> ExpData)
{
	//This function will return a value for the displacement of the AFM 
	//based on experimetnal data.
	int N = ExpData[section].size(); //size of the experimental data
	double* Insert = new double[N];
	double* Extract = new double[N];
	delete[] Insert;
	delete[] Extract;
	return LT_d;
}
//Generates the data for the AFM probe displacement when experimental data is included
double* LT_Set(double* LT, double* time, double L3, double L4, double dt, double speed, double pause, double timeApp, double timeRet, double initDist, int Nt)
{
	//This functions will set function that will be imprinting on the surface. 
	//To begin with assigning for the approach
	double c1;
	double t = 0.0;
	double D = -80.0;
	double val = L3 + L4 - D;
	c1 = timeApp + pause;
	for (int i = 0; i < Nt; i++)
	{
		if (time[i] < timeApp)
		{
			LT[i] = (-speed * time[i]) + initDist;
			LT[i] = L3 + L4 - D + LT[i];
			val = LT[i];
		}
		else if (time[i] > c1)
		{
			LT[i] = (speed * t); 
			LT[i] = val + LT[i];
			t += dt;
		}
		else
		{
			LT[i] = LT[i - 1];
			val = LT[i]; 
		}
	}
	
	return LT;
}
//The rhs of the ODE for the integrin binding in the layer between the cell and AFM probe
double b_rhs(double d, double B, double eta, double** b, double t, int i, int tn)
{
	double rhs_b, c1, c2;
	double h = 10.0; 
	double theta = 2.0;
	double k = 5.0;
	c1 = h * k_bind(d) * (1.0 - B);
	c2 = k_unbind_U(d) * b[i][tn];
	rhs_b = c1 - c2;
	return rhs_b;
}
//Solves the ODE for the inegrin binding between the cell and the AFM probe
double** solve_Diff_upper(double** b, double* x, double* t, double eta, double dx, double dt, double L, int Nx, int tn, int Nt)
{
	double* vec = new double[Nx];
	double d;
	double B = 0.0;
	double h = dt / 720.0;
	double bk1, bk2, bk3, bk4, bP, bF1, bF2, bF3, bF4, bF5, bval, bval_RK;
	for (int i = 0; i < Nx; i++)
	{
		vec[i] = b[i][tn];
	}
	B = trapz(vec, dx, Nx);
	for (int i = 0; i < Nx; i++)
	{
		d = two_dimensional_dist(x[i], L);
		if (tn < 5)
		{
			bval_RK = b[i][tn];
			bk1 = b_rhs(d, B, eta, b, t[tn], i, tn);
			b[i][tn] = b[i][tn] + (dt * (bk1 / 2.0));
			bk2 = b_rhs(d, B, eta, b, t[tn], i, tn);
			b[i][tn] = b[i][tn] + (dt * (bk2 / 2.0));
			bk3 = b_rhs(d, B, eta, b, t[tn], i, tn);
			b[i][tn] = b[i][tn] + (dt * (bk3));
			bk4 = b_rhs(d, B, eta, b, t[tn], i, tn);
			b[i][tn + 1] = bval_RK + ((dt / 6.0) * (bk1 + 2.0 * bk2 + 2.0 * bk3 + bk4));
		}
		else
		{
			bF1 = 1901.0 * b_rhs(d, B, eta, b, t[tn], i, tn);
			bF2 = -2774.0 * b_rhs(d, B, eta, b, t[tn - 1], i, tn - 1);
			bF3 = 2616.0 * b_rhs(d, B, eta, b, t[tn - 2], i, tn - 2);
			bF4 = -1274.0 * b_rhs(d, B, eta, b, t[tn - 3], i, tn - 3);
			bF5 = 251.0 * b_rhs(d, B, eta, b, t[tn - 4], i, tn - 4);
			bP = b[i][tn] + (h * (bF1 + bF2 + bF3 + bF4 + bF5));

			bval = b[i][tn];
			b[i][tn] = bP;

			bF1 = b_rhs(d, B, eta, b, t[tn], i, tn);
			b[i][tn] = bval;
			bF2 = b_rhs(d, B, eta, b, t[tn], i, tn);
			bF3 = b_rhs(d, B, eta, b, t[tn - 1], i, tn - 1);
			bF4 = b_rhs(d, B, eta, b, t[tn - 2], i, tn - 2);
			bF5 = b_rhs(d, B, eta, b, t[tn - 3], i, tn - 3);

			b[i][tn + 1] = b[i][tn] + (h * (251.0 * bF1 + 646.0 * bF2 - 264.0 * bF3 + 106.0 * bF4 - 19.0 * bF5));
		}
	}
	return b;
}
//The right-hs of the ODE for the integrin binding to the biomaterial in the contracted state
double c_rhs(double d, double E_C, double gamma_bar, double theta_bar,  double l1, double l2, double** c, double** e, int lambda1, int lambda2, double t, int i, int tn)
{
	double rhs_c, c1, c2, c3, c4;
	c1 = gamma_bar * k_bind(d) * (1.0 - E_C);
	c2 = theta_bar * k_unbind(d) * c[i][tn];
	c3 = k_e_c(lambda2, l2) * e[i][tn];
	c4 = k_c_e(lambda1, l1) * c[i][tn];
	rhs_c = c1 - c2 + c3 - c4;
	return rhs_c;
}
//The right-hs of the ODE for the integrin binding to the biomaterial in the extended state 
double e_rhs(double d, double E_C, double gamma_hat, double theta_hat, double l1, double l2, double** c, double** e, int lambda1, int lambda2, double t, int i, int tn)
{
	double rhs_e, c1, c2, c3, c4;
	c1 = gamma_hat * k_bind(d) * (1.0 - E_C);
	c2 = theta_hat * k_unbind(d) * e[i][tn];
	c3 = k_e_c(lambda2, l2) * e[i][tn];
	c4 = k_c_e(lambda1, l1) * c[i][tn];
	rhs_e = c1 - c2 - c3 + c4;
	return rhs_e;
}
//Solves the ODE for the integrin layer inbetween the cell and the biomaterial 
double** solve_Diff_lower(double** e, double** c, double* x, double* t, int lambda1, int lambda2, double dx, double dt, double L, double a, double b, double gamma_bar, double gamma_hat, double theta_bar, double theta_hat, double l1, double l2, int Nx, int tn, int Nt)
{
	double* vec = new double[Nx];
	double d;
	double E_C = 0.0;
	double h = dt / 720.0;
	double ck1, ck2, ck3, ck4, cP, cF1, cF2, cF3, cF4, cF5, cval, cval_RK;
	double ek1, ek2, ek3, ek4, eP, eF1, eF2, eF3, eF4, eF5, eval, eval_RK;
	for (int i = 0; i < Nx; i++)
	{
		vec[i] = c[i][tn] + e[i][tn];
	}
	E_C = trapz(vec, dx, Nx);
	for (int i = 0; i < Nx; i++)
	{
		d = two_dimensional_dist(x[i], L);
		if (tn < 5)
		{
			cval_RK = c[i][tn];
			ck1 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			c[i][tn] = c[i][tn] + (dt * (ck1 / 2.0));
			ck2 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			c[i][tn] = c[i][tn] + (dt * (ck2 / 2.0));
			ck3 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			c[i][tn] = c[i][tn] + (dt * (ck3));
			ck4 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			c[i][tn + 1] = cval_RK + ((dt / 6.0) * (ck1 + 2.0 * ck2 + 2.0 * ck3 + ck4));

			eval_RK = e[i][tn];
			ek1 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			e[i][tn] = e[i][tn] + (dt * (ek1 / 2.0));
			ek2 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			e[i][tn] = e[i][tn] + (dt * (ek2 / 2.0));
			ek3 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			e[i][tn] = e[i][tn] + (dt * (ek3));
			ek4 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			e[i][tn + 1] = eval_RK + ((dt / 6.0) * (ek1 + 2.0 * ek2 + 2.0 * ek3 + ek4));
		}
		else
		{
			cF1 = 1901.0 * c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			cF2 = -2774.0 * c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn - 1], i, tn - 1);
			cF3 = 2616.0 * c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn - 2], i, tn - 2);
			cF4 = -1274.0 * c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn - 3], i, tn - 3);
			cF5 = 251.0 * c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn - 4], i, tn - 4);
			cP = c[i][tn] + (h * (cF1 + cF2 + cF3 + cF4 + cF5));

			eF1 = 1901.0 * e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			eF2 = -2774.0 * e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn - 1], i, tn - 1);
			eF3 = 2616.0 * e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn - 2], i, tn - 2);
			eF4 = -1274.0 * e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn - 3], i, tn - 3);
			eF5 = 251.0 * e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn - 4], i, tn - 4);
			eP = e[i][tn] + (h * (eF1 + eF2 + eF3 + eF4 + eF5));

			cval = c[i][tn];
			c[i][tn] = cP;

			eval = e[i][tn];
			e[i][tn] = eP;

			cF1 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			c[i][tn] = cval;
			cF2 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			cF3 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn - 1], i, tn - 1);
			cF4 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn - 2], i, tn - 2);
			cF5 = c_rhs(d, E_C, gamma_bar, theta_bar, l1, l2, c, e, lambda1, lambda2, t[tn - 3], i, tn - 3);

			eF1 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			e[i][tn] = eval;
			eF2 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn], i, tn);
			eF3 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn - 1], i, tn - 1);
			eF4 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn - 2], i, tn - 2);
			eF5 = e_rhs(d, E_C, gamma_hat, theta_hat, l1, l2, c, e, lambda1, lambda2, t[tn - 3], i, tn - 3);

			c[i][tn + 1] = c[i][tn] + (h * (251.0 * cF1 + 646.0 * cF2 - 264.0 * cF3 + 106.0 * cF4 - 19.0 * cF5));
			e[i][tn + 1] = e[i][tn] + (h * (251.0 * eF1 + 646.0 * eF2 - 264.0 * eF3 + 106.0 * eF4 - 19.0 * eF5));
		}
	}
	delete[] vec;
	return c, e;
}
//Calculates the force in the same manner as the experimental data
double Force_calc(double k, double Ln, double Ln_d)
{
	double F;
	F = k * (Ln_d - Ln);
	F = -F;
	return F;
}
