/*Useful functions for debugging code*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "Useful.h"
//Finds the maximum value in a 1D array
double abs_max_val_1d(double* vec, int N)
{
	double max=vec[0];
	for (int i = 1; i < N; i++)
	{
		if (abs(vec[i]) > abs(vec[i - 1]))
		{
			max = abs(vec[i]);
		}
	}
	return max;
}
//Destroys the memory of a 2D array of the type integer
void destroyTwoDimenArray_int(int** ptr, int row, int col)
{
	for (int i = 0; i < row; i++)
	{
		delete[] ptr[i];
	}
	delete[] ptr;
}
//Destroys the memory of a 2D array of the data type double
void destroyTwoDimenArray_double(double** ptr, int row, int col)
{
	for (int i = 0; i < row; i++)
	{
		delete[] ptr[i];
	}
	delete[] ptr;
}
//sgn function 
double sgn(double v)
{
	if (v < 0) return -1.0;
	if (v > 0) return 1.0;
	return 0;
}
//Fills in the time domain for an array with uniform distribution
double* tDomain(double* t, double at, double bt, double dt, int Nt)
{
	t[0] = at;
	for (int i = 0; i < Nt; i++)
	{
		t[i + 1] = t[i] + dt;
	}
	return t;
}
//Fills in the spatial domain for an array with unifrom distribution
double* xDomain(double* x, double ax, double bx, double dx, int Nx)
{
	x[0] = ax;
	for (int i = 0; i < Nx - 1; i++)
	{
		x[i + 1] = x[i] + dx;
	}
	return x;
}
//Integrates a function using the trapezoidal method
double trapz(double* f, double dx, int N)
{
	double I;
	double c1 = (dx / 2.0);
	I = f[0] + f[N - 1];
	for (int i = 1; i < N-1; i++)
	{
		I += 2.0 * f[i];
	}
	I = c1 * I;
	return I;
}
//Calculate a distance using the Pythagorus method
double two_dimensional_dist(double x, double L)
{
	double y, c1, c2, c3;
	c1 = pow(x, 2.0);
	c2 = pow(L, 2.0);
	c3 = c1 + c2;
	y = pow(c3, 0.5);
	return y;
}
//Finds the length of an array
int getLengthFromFile(const std::string& filename) 
{
	std::ifstream file(filename);

	int length = 0;
	std::string line;
	while (std::getline(file, line)) {
		length++;
	}

	file.close();
	return length;
}