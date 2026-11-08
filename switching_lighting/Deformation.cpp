#include <iostream>
#include <cmath>
#include "Deformation.h"
#include "Useful.h"
//See the paper for details: Function calculates a parameter used in the 
//master equation
double gamma(double* k, int N)
{
	double g;
	double dN = N;
	double c1;
	double vec = 0.0;
	for (int i = 0; i < N - 1; i++)
	{
		vec = vec + (1.0 / k[i]);
	}
	c1 = k[N - 1] * vec;
	g = 1.0 + c1;
	return g;
}
//For solving the deformation master equation
double L_master(double* t, double* k, double LT, double LT_d, double L_n, int n, int N, int tn)
{
	double L_prime, c1, c2, c3;
	c1 = k[N - 1] / k[n];
	c2 = (LT_d - LT);
	c3 = (c1 * c2) / gamma(k, N);
	L_prime = L_n + c3;
	return L_prime;
}