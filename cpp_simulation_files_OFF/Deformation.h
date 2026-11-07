#pragma once
//See the paper for details: Function calculates a parameter used in the 
//master equation
double gamma(double* k, int N);
//For solving the deformation master equation
double L_master(double* t, double* k, double LT, double LT_d, double L_n, int n, int N, int tn);